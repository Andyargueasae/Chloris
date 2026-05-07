#!/usr/bin/env python3
"""Generate protein MDX pages from protein JSON metadata.

Reads all files from ``src/proteins_jsons`` and writes generated pages to
``src/content/docs/proteins``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
JSON_DIR = ROOT / "src" / "proteins_jsons"
OUT_DIR = ROOT / "src" / "content" / "docs" / "proteins"
PUBLIC_ALIGNMENTS_DIR = ROOT / "public" / "alignments"
PUBLIC_STRUCTURE_DIR = ROOT / "public" / "structure"


def parse_bracket_list(value: Any) -> list[str]:
	"""Parse values like "[A,B]" into ["A", "B"]."""
	if value is None:
		return []
	if isinstance(value, list):
		return [str(item).strip() for item in value if str(item).strip()]
	text = str(value).strip()
	if not text:
		return []
	if text.startswith("[") and text.endswith("]"):
		text = text[1:-1]
	if not text.strip():
		return []
	items: list[str] = []
	for item in text.split(","):
		cleaned = item.strip().strip('"').strip("'")
		if cleaned:
			items.append(cleaned)
	return items


def cleaned_text(value: Any) -> str:
	if value is None:
		return ""
	text = str(value).strip()
	if not text:
		return ""
	if text.lower() in {"none", "null", "nan"}:
		return ""
	return text


def normalize_alignment_url(protein: str, raw_alignment_path: Any) -> str:
	candidate_name = None
	if raw_alignment_path:
		candidate_name = Path(str(raw_alignment_path)).name
	if not candidate_name:
		candidate_name = f"{protein}.aln"

	preferred = PUBLIC_ALIGNMENTS_DIR / candidate_name
	if preferred.exists():
		return f"../../alignments/{candidate_name}"

	fallback_name = f"{protein}.aln"
	fallback = PUBLIC_ALIGNMENTS_DIR / fallback_name
	if fallback.exists():
		return f"../../alignments/{fallback_name}"

	return f"../../alignments/{candidate_name}"


def normalize_structure_url(protein: str, raw_structure_path: Any) -> str:
	candidate_name = None
	if raw_structure_path:
		candidate_name = Path(str(raw_structure_path)).name
	protein_lower = protein.lower()

	expected_name = f"{protein_lower}_model.cif"
	expected_file = PUBLIC_STRUCTURE_DIR / expected_name
	if expected_file.exists():
		return f"../../structure/{expected_name}"

	if candidate_name and (PUBLIC_STRUCTURE_DIR / candidate_name).exists():
		return f"../../structure/{candidate_name}"

	if candidate_name:
		return f"../../structure/{candidate_name}"

	return f"../../structure/{expected_name}"


def yaml_double_quoted(value: str) -> str:
	escaped = value.replace("\\", "\\\\").replace('"', '\\"')
	return f'"{escaped}"'


def md_inline(value: Any) -> str:
	text = str(value)
	text = text.replace("|", "\\|")
	return text


def build_protein_mdx(protein: str, data: dict[str, Any]) -> str:
	title = cleaned_text(data.get("protein")) or protein
	full_name = cleaned_text(data.get("full name"))
	descriptions = cleaned_text(data.get("descriptions"))
	uniref90 = cleaned_text(data.get("Uniref90 representative"))
	uniprot_link = cleaned_text(data.get("UniProt Representative link"))
	n_seqs = data.get("n_seqs")

	pfam_domains = parse_bracket_list(data.get("pfam domains"))
	pfam_links = parse_bracket_list(data.get("pfam_links"))

	alignment_url = normalize_alignment_url(protein, data.get("alignment_path"))
	structure_url = normalize_structure_url(protein, data.get("structure_path"))

	description_text = full_name or f"{title} alignment and predicted structure"
	lines: list[str] = [
		"---",
		f"title: {yaml_double_quoted(title)}",
		f"description: {yaml_double_quoted(description_text)}",
		f"slug: {yaml_double_quoted(f'proteins/{protein}')}",
		"---",
		"",
		"import MsaBrowser from '../../../components/MSA.astro';",
		"import ProteinStructureViewer from '../../../components/ProteinViewer.astro';",
		"import { LinkButton } from '@astrojs/starlight/components';",
		"",
		"# Protein Information",
		"",
	]

	if full_name:
		lines.extend([f"**Full name:** {full_name}", ""])
	if descriptions:
		lines.extend([descriptions, ""])

	lines.extend(["## Metadata", ""])

	if uniref90:
		lines.append(f"- **UniRef90 representative:** `{uniref90}`")
	if uniprot_link:
		lines.append(f"- **UniProt entry:** [{uniprot_link}]({uniprot_link})")
	if n_seqs is not None:
		try:
			n_seqs_int = int(float(n_seqs))
			lines.append(f"- **Number of aligned sequences:** {n_seqs_int}")
		except (TypeError, ValueError):
			lines.append(f"- **Number of aligned sequences:** {n_seqs}")

	if pfam_domains:
		lines.append("- **Pfam domains:**")
		for idx, domain in enumerate(pfam_domains):
			if idx < len(pfam_links):
				lines.append(f"  - [{domain}]({pfam_links[idx]})")
			else:
				lines.append(f"  - {domain}")
	lines.append("")

	lines.extend(
		[
			"# Multi-Sequence Alignment",
			"",
			f"<MsaBrowser viewerId=\"{protein}Viewer\" fastaUrl=\"{alignment_url}\" colorSchema=\"clustal\" hasConsensus={{true}} exportName=\"{protein}_alignment.fasta\" />",
			"",
			"# Structure Visualization",
			"",
			f"<ProteinStructureViewer cifUrl=\"{structure_url}\" height=\"800px\" />",
			"",
			"# Downloads",
			"",
			f"<LinkButton href=\"{alignment_url}\" icon=\"download\">{protein}.aln</LinkButton>",
			f"<LinkButton href=\"{structure_url}\" icon=\"download\">{protein}.cif</LinkButton>",
			"",
			"# Raw Metadata",
			"",
			"| Field | Value |",
			"| --- | --- |",
		]
	)

	for key, value in data.items():
		if isinstance(value, (dict, list)):
			value_text = json.dumps(value, ensure_ascii=True)
		else:
			value_text = str(value)
		lines.append(f"| {md_inline(key)} | {md_inline(value_text)} |")

	lines.append("")
	return "\n".join(lines)


def write_protein_pages(dry_run: bool = False) -> tuple[int, int]:
	if not JSON_DIR.exists():
		raise FileNotFoundError(f"Input directory does not exist: {JSON_DIR}")

	OUT_DIR.mkdir(parents=True, exist_ok=True)

	generated = 0
	skipped = 0

	for json_path in sorted(JSON_DIR.glob("*.json")):
		protein = json_path.stem
		try:
			with json_path.open("r", encoding="utf-8") as infile:
				data = json.load(infile)
		except json.JSONDecodeError as exc:
			print(f"Skipping {json_path.name}: invalid JSON ({exc})")
			skipped += 1
			continue

		output_path = OUT_DIR / f"{protein}.mdx"
		mdx_content = build_protein_mdx(protein, data)

		if dry_run:
			print(f"Would write {output_path.relative_to(ROOT)}")
		else:
			output_path.write_text(mdx_content, encoding="utf-8")
			print(f"Wrote {output_path.relative_to(ROOT)}")
		generated += 1

	return generated, skipped


def main() -> None:
	parser = argparse.ArgumentParser(description="Build protein MDX pages from protein JSON files.")
	parser.add_argument("--dry-run", action="store_true", help="List files without writing output.")
	args = parser.parse_args()

	generated, skipped = write_protein_pages(dry_run=args.dry_run)
	print(f"Done. Generated {generated} page(s), skipped {skipped} file(s).")


if __name__ == "__main__":
	main()
