# PRISCREEN panel format (for PRISTINE)

PRISCREEN identifies reads of PRISTINE amplicons against a *panel*: one FASTA
record per amplicon variant, labelled with the species it comes from. PRISCREEN
can build the panel from `<group>/primers/amplicons/<locus>_amplicons.fasta`
(`scripts/pristine_to_panel.py`), but a panel written by PRISTINE itself is
used directly and is preferred, because PRISTINE knows the species of every
genome and which genomes are targets.

## File

- `<output>/<group>/primers/priscreen_panel.fasta`: the panel of one group.
- `<output>/priscreen_panel.fasta`: all groups of the run merged (written at
  the end of the run from the group files that exist).

PRISCREEN looks for the file at an output root, in a group directory
(`primers/priscreen_panel.fasta`) or in a `primers/` directory, and falls back
to the `amplicons/*_amplicons.fasta` files of any group without one.

## Records

```text
>locus|species=Genus_species|genome=Genome_id|role=target
ACGT...
```

- `locus`: the locus name PRISTINE uses for the primer pair
  (`<contig>_<start>_<end>` on the reference), without `|`. In the merged
  file a locus name that occurs in two groups is prefixed with the group name
  (`<group>_<locus>`) in both.
- `species=`: species of the genome, as PRISTINE defines it (file name without
  the trailing `_<i>`; for a multi-species target group the species of the
  genome, not the group name). No spaces: `Yersinia_pestis`.
- `genome=`: the genome id (file stem). Optional for PRISCREEN, useful for
  tracing a hit.
- `role=`: `target` or `non-target` (`background` for an amplifying background
  product, if PRISTINE wants to include those). Optional; PRISCREEN ignores
  unknown `key=value` fields.
- Sequence: the predicted amplicon in the orientation of the reference locus,
  upper case, one record per *distinct* sequence of a (locus, species): two
  genomes of a species with the identical amplicon give one record (keep the
  first genome id).
- Fields are separated by `|`; `locus` and `species=` are mandatory; the header
  must not contain spaces before the last field (PRISCREEN uses the text up to
  the first whitespace).

Example from the toy data of PRISCREEN (`data/reference_amplicons.fasta`
uses the minimal form):

```text
>amp1|species=SpeciesA
>amp1|species=SpeciesB
```
