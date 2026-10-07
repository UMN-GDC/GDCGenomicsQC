# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- NMIND Bronze compliance: Zenodo metadata (`.zenodo.json`), GitHub Actions CI workflow
- SNP heritability estimation for related populations: fastGWA (GCTA) and PredLMM-ACE (SAONLIB)
- New `snpHeritRelated` config block for unified related/unrelated heritability runs
- Wrapper scripts: `run_predlmm_ace.py`, `parse_fastgwa.py`
- Conda environment `predlmmAce.yml` with gcta, predlmm-ace, and dependencies
- Container images for all environments via GHCR

### Changed
- Updated README with CI, DOI, and release badges
- Extended config schema with `snpHeritRelated` parameters
- Documentation index with badges

### Fixed
- Internal PCA method selection for related vs unrelated GRM

## [1.0.0] - 2026-10-07

### Added
- Initial release of GDCGenomicsQC pipeline
- Snakemake-based modular QC workflow
- Relatedness estimation (KING, PRIMUS)
- Global ancestry classification (PCA, UMAP, RFMix local ancestry)
- SNP heritability estimation (AdjHE, GCTA, PredLMM, SWD, COMBAT, Covbat)
- Phenotype simulation with controlled heritability
- Polygenic risk scoring (single and multi-ancestry methods: C+T, PRSice, PRS-CS, LDpred2, lassosum2, PRS-CSx, PROSPER)
- Containerized execution via Apptainer/Singularity
- HPC profiles (MSI, sandbox, interactive)
- Automated report generation with Sphinx/ReadTheDocs
- Comprehensive tutorial documentation

### Changed
- N/A (initial release)

### Deprecated
- N/A

### Removed
- N/A

### Fixed
- N/A

### Security
- N/A

---

## Release Process

1. Update version in relevant files (none currently, using git tags)
2. Update this CHANGELOG.md with release notes
3. Create git tag: `git tag -a v1.0.0 -m "Release v1.0.0"`
4. Push tag: `git push origin v1.0.0`
5. GitHub Actions will:
   - Run full CI pipeline
   - Build and push container images to GHCR
   - Create GitHub Release with auto-generated changelog
   - Deploy documentation to GitHub Pages
6. Manually create Zenodo record from GitHub release (or enable auto-sync)

---

## Version History

| Version | Date | Notes |
|---------|------|-------|
| 1.0.0 | 2026-10-07 | Initial public release |