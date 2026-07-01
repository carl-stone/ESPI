import { spawnSync } from "node:child_process";
import path from "node:path";
import { fileURLToPath } from "node:url";
import type { ExtensionAPI } from "@oh-my-pi/pi-coding-agent";

type ToolEvent = {
  toolName?: unknown;
  input?: unknown;
};

type ToolCallDecision = {
  block: true;
  reason: string;
};

const adapterDir = path.dirname(fileURLToPath(import.meta.url));
const policyPath = path.resolve(adapterDir, "../policies/generated-files.sh");

function isRecord(input: unknown): input is Record<string, unknown> {
  return Boolean(input && typeof input === "object" && !Array.isArray(input));
}

const WRITE_TOOLS: Record<string, true> = {
  apply_patch: true,
  ast_edit: true,
  edit: true,
  write: true,
};

function directPaths(input: unknown): string[] {
  if (!isRecord(input)) return [];

  const paths: string[] = [];
  for (const key of ["path", "file_path", "filePath", "filename"]) {
    const value = input[key];
    if (Array.isArray(value)) {
      for (const item of value) {
        if (typeof item === "string" && item.length > 0) paths.push(item);
      }
    } else if (typeof value === "string" && value.length > 0) {
      paths.push(value);
    }
  }

  return paths;
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
  if (!text) return [];

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

function writeRiskPaths(event: ToolEvent): string[] {
  const toolName = String(event.toolName ?? "").toLowerCase();
  const input = event.input;

  if (!WRITE_TOOLS[toolName]) return [];

  const seen = new Set<string>();
  const paths: string[] = [];
  for (const candidate of [...directPaths(input), ...patchPaths(textPayload(input))]) {
    if (candidate.length === 0 || seen.has(candidate)) continue;
    seen.add(candidate);
    paths.push(candidate);
  }

  return paths;
}

function runPolicy(paths: string[], cwd: string): ToolCallDecision | undefined {
  if (paths.length === 0) return undefined;

  const args = ["--operation", "write", ...paths.flatMap((item) => ["--path", item])];
  const result = spawnSync(policyPath, args, {
    cwd,
    encoding: "utf8",
  });

  const output = [result.stderr, result.stdout]
    .filter((item): item is string => Boolean(item && item.trim()))
    .join("\n")
    .trim();

  if (result.status === 0) return undefined;

  if (result.status === 2) {
    return {
      block: true,
      reason: output || "Generated-file policy blocked this edit.",
    };
  }

  return {
    block: true,
    reason: output
      ? `Generated-file policy failed; blocking as a precaution. ${output}`
      : "Generated-file policy failed; blocking as a precaution.",
  };
}

export default function generatedFileGuard(pi: ExtensionAPI): void {
  pi.setLabel("Generated File Guard");

  pi.on("tool_call", async (event, ctx) => {
    const paths = writeRiskPaths(event as ToolEvent);
    return runPolicy(paths, ctx.cwd);
  });
}
