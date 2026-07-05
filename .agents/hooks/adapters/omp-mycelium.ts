import { spawnSync } from "node:child_process";
import fs from "node:fs";
import path from "node:path";
import type { ExtensionAPI } from "@oh-my-pi/pi-coding-agent";

type ToolEvent = {
  toolName?: unknown;
  input?: unknown;
  content?: unknown;
  isError?: unknown;
};

type SessionDecision = {
  decision: "block";
  reason: string;
};

type HookJson = {
  decision?: unknown;
  reason?: unknown;
  additionalContext?: unknown;
};

type ToolResultDecision = {
  content?: unknown[];
  details?: unknown;
  isError?: boolean;
};

const WRITE_TOOLS: Record<string, true> = {
  apply_patch: true,
  ast_edit: true,
  edit: true,
  write: true,
};

function isRecord(input: unknown): input is Record<string, unknown> {
  return Boolean(input && typeof input === "object" && !Array.isArray(input));
}

function ensureDir(filePath: string): void {
  fs.mkdirSync(path.dirname(filePath), { recursive: true });
}

function appendUniqueLine(filePath: string, line: string): void {
  ensureDir(filePath);
  const existing = fs.existsSync(filePath) ? fs.readFileSync(filePath, "utf8") : "";
  if (existing.split(/\r?\n/).includes(line)) return;
  fs.appendFileSync(filePath, `${line}\n`);
}

function appendLine(filePath: string, line: string): void {
  ensureDir(filePath);
  fs.appendFileSync(filePath, `${line}\n`);
}

function repoRoot(cwd: string): string {
  const result = spawnSync("git", ["rev-parse", "--show-toplevel"], {
    cwd,
    encoding: "utf8",
  });
  if (result.status !== 0) return cwd;
  const root = result.stdout.trim();
  return root.length > 0 ? root : cwd;
}

function runHook(cwd: string, hookPath: string, payload: unknown): HookJson | undefined {
  if (!fs.existsSync(path.resolve(cwd, hookPath))) return undefined;

  const result = spawnSync(hookPath, {
    cwd,
    input: JSON.stringify(payload),
    encoding: "utf8",
    shell: false,
  });

  const text = result.stdout.trim();
  if (text.length === 0) return undefined;

  try {
    return JSON.parse(text) as HookJson;
  } catch {
    return undefined;
  }
}

function hookAdditionalContext(result: HookJson | undefined): string | undefined {
  if (!result || typeof result.additionalContext !== "string") return undefined;
  const text = result.additionalContext.trim();
  return text.length > 0 ? text : undefined;
}

function appendHookContext(event: ToolEvent, result: HookJson | undefined): ToolResultDecision | undefined {
  const text = hookAdditionalContext(result);
  if (!text) return undefined;

  const existing = Array.isArray(event.content) ? event.content : [];
  return {
    content: [...existing, { type: "text", text }],
  };
}

function runSync(cwd: string): void {
  const script = path.resolve(cwd, "tools/sync-mycelium-skills-core.py");
  if (!fs.existsSync(script)) return;
  spawnSync("python3", [script, "--quiet"], { cwd, encoding: "utf8" });
}

function stringValues(input: unknown, keys: string[]): string[] {
  if (!isRecord(input)) return [];

  const values: string[] = [];
  for (const key of keys) {
    const value = input[key];
    if (Array.isArray(value)) {
      for (const item of value) {
        if (typeof item === "string" && item.length > 0) values.push(item);
      }
    } else if (typeof value === "string" && value.length > 0) {
      values.push(value);
    }
  }
  return values;
}

function directPaths(input: unknown): string[] {
  return stringValues(input, ["path", "file_path", "filePath", "filename"]);
}

function textPayload(input: unknown): string {
  if (typeof input === "string") return input;
  if (!isRecord(input)) return "";

  for (const key of ["patch", "command", "input", "content", "text"]) {
    const value = input[key];
    if (typeof value === "string" && value.length > 0) return value;
  }
  return "";
}

function patchPaths(text: string): string[] {
  if (text.length === 0) return [];

  const paths: string[] = [];
  const add = (candidate: string | undefined) => {
    if (candidate && candidate !== "/dev/null") paths.push(candidate.trim());
  };

  for (const line of text.split(/\r?\n/)) {
    let match = line.match(/^\*\*\* (?:Add|Update|Delete) File: (.+)$/);
    if (match) {
      add(match[1]);
      continue;
    }

    match = line.match(/^\*\*\* Rename (?:from|to): (.+)$/);
    if (match) {
      add(match[1]);
      continue;
    }

    match = line.match(/^\[([^#\]\r\n]+)#[0-9A-Fa-f]{4}\]$/);
    if (match) {
      add(match[1]);
      continue;
    }

    match = line.match(/^--- a\/(.+)$/);
    if (match) {
      add(match[1]);
      continue;
    }

    match = line.match(/^\+\+\+ b\/(.+)$/);
    if (match) add(match[1]);
  }

  return paths;
}

function normalizeFilePath(candidate: string, root: string): string {
  if (candidate.startsWith("local://") || candidate.startsWith("artifact://")) return candidate;
  const stripped = candidate.replace(/:(?:raw|\d+(?:[-+]\d+)?(?:,\d+(?:[-+]\d+)?)*)$/, "");
  return path.isAbsolute(stripped) ? stripped : path.resolve(root, stripped);
}

function livingLogPath(candidate: string, root: string): string | undefined {
  const normalized = normalizeFilePath(candidate, root);
  const marker = `${path.sep}.living${path.sep}`;
  const markerIndex = normalized.indexOf(marker);
  if (markerIndex >= 0) return `.living/${normalized.slice(markerIndex + marker.length)}`;
  if (candidate === ".living" || candidate.startsWith(".living/")) return candidate;
  return undefined;
}

function skipActivityPath(candidate: string, root: string): boolean {
  const normalized = normalizeFilePath(candidate, root);
  return normalized.includes(`${path.sep}.living${path.sep}`) || normalized.includes(`${path.sep}.claude${path.sep}`);
}

function isMyceliumMaintenanceCommand(command: string): boolean {
  const normalized = command.replace(/\\/g, "/").trim();
  if (/[;&|]/.test(normalized)) return false;

  const pythonRunner = String.raw`(?:python3?|uv\s+run\s+python3?)`;
  const optionalPrefix = String.raw`(?:[^\s]*/)?`;
  const coreScripts = [
    "generate_index",
    "validate_structure",
    "recall_lessons",
    "detect_recurrence",
    "crystallize_findings",
    "init_knowledge",
    "migrate_existing_repos",
    "install_convention",
    "init_repo",
  ].join("|");
  const script = String.raw`(?:${optionalPrefix}skills/core/scripts/(?:${coreScripts})\.py|${optionalPrefix}tools/sync-mycelium-skills-core\.py)`;
  const pattern = new RegExp(String.raw`^${pythonRunner}\s+${script}(?:\s+.*)?$`);
  return pattern.test(normalized);
}

function sessionId(input: unknown): string | undefined {
  if (!isRecord(input)) return undefined;
  for (const key of ["session_id", "sessionId"]) {
    const value = input[key];
    if (typeof value === "string" && value.length > 0) return value;
  }
  return undefined;
}


function writeReadAccess(toolPath: string, root: string): void {
  const relPath = livingLogPath(toolPath, root);
  if (!relPath) return;

  const timestamp = new Date().toISOString().slice(0, 19);
  appendLine(path.join(root, ".claude/mycelium-read-access.log"), `${timestamp} ${relPath}`);
}

function writeActivity(paths: string[], root: string): void {
  const activityPath = path.join(root, ".claude/mycelium-session-activity.tmp");
  const usefulPaths = paths
    .map((item) => normalizeFilePath(item, root))
    .filter((item) => !skipActivityPath(item, root));

  for (const filePath of usefulPaths) appendUniqueLine(activityPath, filePath);
  if (usefulPaths.length === 0) return;

  const reminderPath = path.join(root, ".claude/mycelium-reminded.tmp");
  if (!fs.existsSync(reminderPath)) {
    ensureDir(reminderPath);
    fs.writeFileSync(reminderPath, `${Math.floor(Date.now() / 1000)}\n`);
  }
}

function markdownSection(markdown: string, heading: string): string {
  const marker = `## ${heading}\n`;
  const start = markdown.indexOf(marker);
  if (start < 0) return "";
  const sectionStart = start + marker.length;
  const rest = markdown.slice(sectionStart);
  const nextHeading = rest.search(/^## /m);
  const section = nextHeading >= 0 ? rest.slice(0, nextHeading) : rest;
  return section.trim();
}

function writeFallbackLastSession(root: string): void {
  const claudeDir = path.join(root, ".claude");
  const activityPath = path.join(claudeDir, "mycelium-session-activity.tmp");
  const lastSessionPath = path.join(claudeDir, "last-session.md");
  const existing = fs.existsSync(lastSessionPath) ? fs.readFileSync(lastSessionPath, "utf8") : "";
  const requiredSections = ["What was worked on", "Key decisions made", "Blockers & surprises", "Current state", "Next steps"];
  const hasFiveSections = requiredSections.every((heading) => existing.includes(`## ${heading}`));
  if (hasFiveSections) return;

  ensureDir(lastSessionPath);

  const files = fs.existsSync(activityPath)
    ? fs.readFileSync(activityPath, "utf8").split(/\r?\n/).filter(Boolean).slice(0, 10)
    : [];
  const branch = spawnSync("git", ["rev-parse", "--abbrev-ref", "HEAD"], {
    cwd: root,
    encoding: "utf8",
  }).stdout.trim() || "unknown";
  const statusLines = spawnSync("git", ["status", "--porcelain"], {
    cwd: root,
    encoding: "utf8",
  }).stdout.split(/\r?\n/).filter(Boolean).length;
  const existingWork = markdownSection(existing, "What was worked on");
  const existingState = markdownSection(existing, "Current state");
  const changed =
    files.length > 0
      ? files.map((file) => `- Modified \`${path.relative(root, file)}\``).join("\n")
      : existingWork || "- No Edit/Write activity was tracked.";
  const uncommitted = statusLines === 1 ? "1 uncommitted change" : `${statusLines} uncommitted changes`;
  const currentState = existingState || `- Branch: \`${branch}\`; ${uncommitted}.`;

  fs.writeFileSync(
    lastSessionPath,
    `# Session resume\n\n## What was worked on\n${changed}\n\n## Key decisions made\n- See \`.living/decisions.md\` for decisions logged during the session.\n\n## Blockers & surprises\n- No blocker was recorded by the OMP Mycelium hook.\n\n## Current state\n${currentState}\n\n## Next steps\n- Continue from the current user request and rerun targeted validation after further changes.\n`,
  );
}

function toolName(event: ToolEvent): string {
  return String(event.toolName ?? "").toLowerCase();
}

export default function myceliumHooks(pi: ExtensionAPI): void {
  pi.setLabel("Mycelium Hooks");

  pi.on("session_start", async (_event, ctx) => {
    const root = repoRoot(ctx.cwd);
    runSync(root);
    runHook(root, ".agents/hooks/adapters/mycelium-health-wrapper.sh", { cwd: root, source: "omp" });
  });

  pi.on("tool_result", async (event, ctx) => {
    const root = repoRoot(ctx.cwd);
    const currentTool = toolName(event as ToolEvent);
    const input = (event as ToolEvent).input;

    if (currentTool === "read") {
      for (const item of directPaths(input)) writeReadAccess(item, root);
      return;
    }

    if (WRITE_TOOLS[currentTool]) {
      writeActivity([...directPaths(input), ...patchPaths(textPayload(input))], root);
      return;
    }

    if (currentTool === "bash") {
      const command = stringValues(input, ["command"])[0];
      if (command && !isMyceliumMaintenanceCommand(command)) {
        const payload = { tool_input: { command }, cwd: ctx.cwd, source: "omp" };
        const result = runHook(root, ".agents/hooks/adapters/mycelium-post-action-wrapper.sh", payload);
        runHook(root, ".agents/hooks/adapters/mycelium-data-tracker-wrapper.sh", payload);
        return appendHookContext(event as ToolEvent, result);
      }
    }
  });

  pi.on("session_stop", async (event, ctx): Promise<SessionDecision | undefined> => {
    const root = repoRoot(ctx.cwd);
    const payload = { cwd: root, source: "omp", session_id: sessionId(event) };
    const result = runHook(root, ".agents/hooks/adapters/mycelium-stop-wrapper.sh", payload);
    runHook(root, "skills/core/hooks/mycelium-data-lineage-stop.sh", payload);
    writeFallbackLastSession(root);

    if (result?.decision === "block" && typeof result.reason === "string") {
      return { decision: "block", reason: result.reason };
    }
    return undefined;
  });
}
