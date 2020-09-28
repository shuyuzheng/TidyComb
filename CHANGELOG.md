## [0.0.8] - 2020.09.28

### Changed

- CalculateMat: Change S score in summary table into three different versions: S_max, S_sum, and S_mean
- CalculateMat/CalculateTemplate/ParCalculateTemplate: Change the table name of "synergy" in output object into "response" to match the table name required by DrugComb

## [0.0.5] - 2019.10.14

### Added

- New function "RIConfidenceInterval" function to calculate confidence interval for RI score

### Changed

- Enable "CalculateSens" function to output "pred" table which contains predicted response, standard deviation from fitted model.
- Fix bug in "PredictResponse", "GetPubPhase" and "smooting"