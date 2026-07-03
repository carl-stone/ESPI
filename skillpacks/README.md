# Skill Packs

External skill repositories cloned here for use by the `skill-bridge` convention pack. These are **inert reference libraries** — never installed as agent skill packs.

## Installed Locally

| Repository | Status | Commit / note |
|------------|--------|---------------|
| `scientific-agent-skills` | cloned | `64d5d22d7d79db2e81bbd5a57f2557fa540ea21c` |
| `bioSkills` | cloned | `1e024ea8547ada12039edbe8197aaa959d97763f` |
| `Autonomous-Science` | not cloned | `https://github.com/arjunrajlaboratory/Autonomous-Science.git` returned `Repository not found`; GitHub search found no matching public repo. |

## Setup

```bash
cd skillpacks/
git clone https://github.com/K-Dense-AI/scientific-agent-skills.git
git clone https://github.com/GPTomics/bioSkills.git
# Await corrected URL before cloning Autonomous-Science.
```

## Updating

```bash
cd skillpacks/scientific-agent-skills && git pull
cd ../bioSkills && git pull
```

## How These Are Used

The `skill-bridge` convention pack (in `.living/conventions/skill-bridge/` or `network/conventions/skill-bridge/`) routes analysis workflows to specific SKILL.md files within these repos. The agent reads one file at a time (~150-200 lines per analysis step), never loading the full repos into context.
