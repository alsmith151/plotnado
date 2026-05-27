# Changelog

## [0.3.1](https://github.com/alsmith151/plotnado/compare/plotnado-v0.3.0...plotnado-v0.3.1) (2026-05-27)


### Features

* Add adjust_plot method to BedSimple for improved axis handling ([f744510](https://github.com/alsmith151/plotnado/commit/f744510da5892d30cf832d3a6ce1eb4fa6c4df51))
* Add autoscaling and autocoloring support to Figure class ([54545d1](https://github.com/alsmith151/plotnado/commit/54545d169fde9f8bb3eedd5fe2400e99f38ed4c8))
* Add basic plotting and simple overlays example notebooks ([7ca75d4](https://github.com/alsmith151/plotnado/commit/7ca75d44d6f7a12e9fd763510c58fb34d111291a))
* Add build_docs.yml for building documentation with GitHub Actions ([e3efdcd](https://github.com/alsmith151/plotnado/commit/e3efdcdbcbf9859346eb85785f13035d028122d8))
* Add Docker definition for building and installing Plotnado with micromamba ([c81eb65](https://github.com/alsmith151/plotnado/commit/c81eb65e8928580bb4b059b2327e8fdecf8bad57))
* add gene filtering functionality to Genes track and update documentation ([6bee424](https://github.com/alsmith151/plotnado/commit/6bee424c12412ba8fc6d60e04b19cdb9476b1911))
* Add iced to dependencies in Docker definition for Plotnado ([1315077](https://github.com/alsmith151/plotnado/commit/13150775eb63674ff1b63bf92cfa034ff328c5b1))
* Add loguru to requirements.txt for improved logging capabilities ([cfed464](https://github.com/alsmith151/plotnado/commit/cfed4645937862af872d8c8987be57ad09cf796c))
* add overlays notebook ([8033d7a](https://github.com/alsmith151/plotnado/commit/8033d7ab51b0f2e2f5c605d5c40959671e063cf4))
* add permissions for issues in release workflow ([f7868eb](https://github.com/alsmith151/plotnado/commit/f7868eb6af0f2b63420b6f26e812682c62a97723))
* Add PlotNado documentation ([d08b9f8](https://github.com/alsmith151/plotnado/commit/d08b9f88187ebe7f5fcdebb59d2960fae5b6db9c))
* Add PlotNado documentation ([8403280](https://github.com/alsmith151/plotnado/commit/8403280f0f3f9cb51b64b7c8dad96f72b034fc07))
* Add pybedtools to requirements.txt for enhanced genomic data processing ([1b903b8](https://github.com/alsmith151/plotnado/commit/1b903b8b72648334ef2d30d9f68a8ff623b2f34b))
* Add pyranges to requirements.txt for enhanced genomic data handling ([b3edb07](https://github.com/alsmith151/plotnado/commit/b3edb07d73409d9cf94faa53998ef24b94e29be4))
* Add support for autocoloring of genes in Figure class ([d5013f3](https://github.com/alsmith151/plotnado/commit/d5013f3dac23df2111e38cc21db0f968221aa9f2))
* added genes track ([94c891a](https://github.com/alsmith151/plotnado/commit/94c891a6813b1b588095628b764a0e9aeab51a69))
* added hg38 genes ([9b0211d](https://github.com/alsmith151/plotnado/commit/9b0211d2d655142cc53f35d55b035054c7204bcc))
* enhance autoscaling functionality and improve overlay track limits handling ([800e26e](https://github.com/alsmith151/plotnado/commit/800e26e4afeddf34e13ffb62f62bb3d792e5c1e0))
* enhance gene label handling and update BigWig diff input types ([edb6301](https://github.com/alsmith151/plotnado/commit/edb630128d41d5f81dc482322bb659b11f5999d2))
* Enhance interval fetching and improve plot labeling in BedSimple and PlotGenes classes ([91c9d25](https://github.com/alsmith151/plotnado/commit/91c9d25e2b1769b69d4531d51848d2a02fd8abc4))
* Enhance interval fetching and improve plot labelling ([421c89f](https://github.com/alsmith151/plotnado/commit/421c89f1198df3fc1e9d4dd395b6b69ab708cbed))
* **Figure:** added autocolor feature ([f733714](https://github.com/alsmith151/plotnado/commit/f733714efca814e096f0422e9e148f2af903ff28))
* Introduce BedBasePlotnado class to enhance file handling ([81ddf06](https://github.com/alsmith151/plotnado/commit/81ddf06c9069faa733d5d601aab5a94e96fd8500))
* Introduce BedBasePlotnado class to enhance file handling in Genes class ([256df8b](https://github.com/alsmith151/plotnado/commit/256df8b7716990194654f39e9e051e2122bba512))
* Update PlotNado installation instructions in index.md ([f3380c8](https://github.com/alsmith151/plotnado/commit/f3380c898f6bff59e9f9ceef0ab58b8b95351b59))
* Update plotnado.cli script path in pyproject.toml ([aa93f52](https://github.com/alsmith151/plotnado/commit/aa93f52d2fdf148950a557a826970457c79d3135))
* Update y_label calculation in genes.py ([3fcb7b3](https://github.com/alsmith151/plotnado/commit/3fcb7b39dc27a5fa7822ddf53f46dbe1fc399191))
* Update y_label calculation in genes.py ([266928f](https://github.com/alsmith151/plotnado/commit/266928f32686f12a384890623e2a5f3b7e18f1af))


### Bug Fixes

* Clean up whitespace and improve error handling in BedSimple ([3cca4df](https://github.com/alsmith151/plotnado/commit/3cca4df25facfd5d009c2207ebb40ce9995c0225))
* Clean up whitespace and improve error handling in BedSimple fetch_data method ([e6b24c2](https://github.com/alsmith151/plotnado/commit/e6b24c27525c9dee0f05919efe5c64089d1345d7))
* Correctly return intervals from fetch_data method in BedSimple ([867d911](https://github.com/alsmith151/plotnado/commit/867d91143e299d9b82518621dfa50280d2fbffa4))
* Correctly return intervals from fetch_data method in BedSimple class ([6a26471](https://github.com/alsmith151/plotnado/commit/6a26471a1e9002b1a3fc068e73bc513205fc9a7b))
* Import iced conditionally for ICE normalization method in MatrixCapcruncher ([bad6836](https://github.com/alsmith151/plotnado/commit/bad683612bbd35fa77d3319beb5ff65822d34873))
* Import iced conditionally for normalization method in MatrixCapcruncher ([0fa542e](https://github.com/alsmith151/plotnado/commit/0fa542ec25b0006808a37dc24f92705b9243c50f))
* Update Genes class to inherit from BedBase instead of BedBasePlo… ([a65ab38](https://github.com/alsmith151/plotnado/commit/a65ab38ed3b73e0e7ca52d7d6abfadc63a1275b1))
* Update Genes class to inherit from BedBase instead of BedBasePlotnado ([296b3b7](https://github.com/alsmith151/plotnado/commit/296b3b7be93bb7021f006f4a04ef3d94e281ad3a))


### Documentation

* added basic example notebook ([0bb719c](https://github.com/alsmith151/plotnado/commit/0bb719c256fd7fc47e4bba3f2c791da1fa7c73a1))
* clean documentation and ensure examples render ([#32](https://github.com/alsmith151/plotnado/issues/32)) ([cce6fbc](https://github.com/alsmith151/plotnado/commit/cce6fbcabfa049812570b66999abce47b96879a6))
* updated documentation ([0aedbc6](https://github.com/alsmith151/plotnado/commit/0aedbc65d4c0e14b370fc24cb69d63d05941682d))
