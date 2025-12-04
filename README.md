
---

# Halo Gas and Stellar Fraction Datasets

This repository contains a collection of halo gas fractions, and in some cases, stellar fractions, compiled from the literature. These datasets are primarily derived from X-ray observations, with total halo masses estimated under the assumption of thermal equilibrium.

> **Disclaimer:** The tables and data are provided *as-is*, without any warranty. Users are responsible for verifying and interpreting the data appropriately. If you use any of these datasets in your research, please cite the original publication.

---

## Summary of Datasets

| Dataset               | Type          | Redshift Range | Halo Mass Range | Reference                                                                                                              |
| --------------------- | ------------- | -------------- | --------------- | ---------------------------------------------------------------------------------------------------------------------- |
| Popesso et al. 2024   | Gas + Stellar | ~0.01–0.5      | 10¹³–10¹⁵ M☉    | [arXiv:2411.16555](https://arxiv.org/pdf/2411.16555)                                                                   |
| Gonzalez et al. 2013  | Stellar       | 0.02–0.3       | 10¹³–10¹⁵ M☉    | [arXiv:1309.3565](https://arxiv.org/pdf/1309.3565.pdf)                                                                 |
| Vikhlinin et al. 2006 | Gas           | 0.0–0.3        | 10¹³–10¹⁵ M☉    | [Journal](https://iopscience.iop.org/article/10.1086/500288/pdf) / [arXiv](https://arxiv.org/pdf/astro-ph/0507092.pdf) |
| Giodini et al. 2009   | Gas + Stellar | 0.1–1          | 10¹³–10¹⁵ M☉    | [arXiv:0904.0448](https://arxiv.org/pdf/0904.0448.pdf)                                                                 |
| Arnaud et al. 2017    | Gas           | 0.0–0.5        | 10¹³–10¹⁵ M☉    | [A&A PDF](https://www.aanda.org/articles/aa/pdf/2007/42/aa8541-07.pdf)                                                 |
| Sun et al. 2008       | Gas           | 0.01–0.1       | 10¹³–10¹⁴ M☉    | [arXiv:0805.2320](https://arxiv.org/pdf/0805.2320.pdf)                                                                 |
| Wicker et al. 2023    | Gas           | 0.0–0.5        | 10¹³–10¹⁵ M☉    | [arXiv:2204.12823](https://arxiv.org/pdf/2204.12823)                                                                   |
| Ettori et al. 2011    | Gas           | 0.0–1.3        | 10¹³–10¹⁵ M☉    | [arXiv:1009.3266](https://arxiv.org/pdf/1009.3266.pdf)                                                                 |
| Kravtsov et al. 2018  | Stellar       | 0.0–1          | 10¹²–10¹⁵ M☉    | [arXiv:1401.7329](https://arxiv.org/pdf/1401.7329.pdf)                                                                 |
| Landry et al. 2013    | Gas           | 0.2–1          | 10¹³–10¹⁵ M☉    | [arXiv:1211.4626](https://arxiv.org/pdf/1211.4626.pdf)                                                                 |
| Reyes et al. 2012     | Stellar       | 0.02–0.3       | 10¹³–10¹⁵ M☉    | [arXiv:1110.4107](https://arxiv.org/pdf/1110.4107.pdf)                                                                 |
| Sanderson et al. 2013 | Gas           | 0.02–0.2       | 10¹³–10¹⁵ M☉    | [arXiv:1212.1613](https://arxiv.org/pdf/1212.1613.pdf)                                                                 |

> Note: The redshift and halo mass ranges are approximate, based on the values reported in the papers.

---

## Usage

Each dataset is provided in a machine-readable format (e.g., CSV, FITS) in the repository. Typical columns include:

* Halo mass estimates
* Gas fraction
* Stellar fraction (if available)
* Associated uncertainties

An example of how to load and plot the datasets is provided in **`plot_gas_fractions.py`**.

---

## Citation

Please cite the original papers corresponding to the dataset(s) you use. Proper attribution ensures credit for the original authors and helps reproducibility in scientific research.

---

## License / Terms of Use

* The datasets in this repository are provided **as-is**, for research and educational purposes only.
* No warranty is given regarding accuracy, completeness, or suitability for any purpose.
* Users are responsible for verifying and interpreting the data.
* When using any of these datasets, please cite the original publication(s).
* Redistribution is allowed, provided this notice and citations are included.

---
