# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

<!-- insertion marker -->
## [v1.6.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.6.0) - 2026-08-31

<small>[Compare with v1.5.0](https://github.com/MolarVerse/PQAnalysis/compare/v1.5.0...v1.6.0)</small>

### Docs

- document analysis performance validation ([f3fee7c](https://github.com/MolarVerse/PQAnalysis/commit/f3fee7ca23e7645e8b3290a1faacb4ad4ce37686) by Josef M. Gallmetzer).

### Bug Fixes

- parse paths and unbracketed lists in PQ input files ([3d9a1ac](https://github.com/MolarVerse/PQAnalysis/commit/3d9a1acdca1a8b577a7a8e7a7280ab508a16a4ae) by Josef M. Gallmetzer).
- validate cell box lengths and angles on construction ([a6b45a1](https://github.com/MolarVerse/PQAnalysis/commit/a6b45a10e195073854e9b0b420f64a51c9f9486a) by Josef M. Gallmetzer).
- make --progress flag match its help text ([2eefbb5](https://github.com/MolarVerse/PQAnalysis/commit/2eefbb5281a0bdc2d7763d5de4ef81b38c70300e) by Josef M. Gallmetzer).
- parse xyz2gen --periodic values from the command line ([e346ac3](https://github.com/MolarVerse/PQAnalysis/commit/e346ac33010dfa90853fa8720f5c75b67a136c7d) by Josef M. Gallmetzer).
- accept --n-molecules flag in add_molecules ([514ed60](https://github.com/MolarVerse/PQAnalysis/commit/514ed606ce958cb426bbcccfc7b5ff6027f60bcc) by Josef M. Gallmetzer).
- keep constructor topology when read() gets none ([b0524e9](https://github.com/MolarVerse/PQAnalysis/commit/b0524e9779ee2fa809178c3700c8c1cefd821891) by Josef M. Gallmetzer).
- raise logged errors independently of the logging level ([b94d2d2](https://github.com/MolarVerse/PQAnalysis/commit/b94d2d2e1b6f0c4068aa5e7a46aaa95474debce5) by Josef M. Gallmetzer).
- accept trajectory files on the traj2qmcfc command line ([6bf8369](https://github.com/MolarVerse/PQAnalysis/commit/6bf83692430cc75f99849517df5c3007e746a62a) by Josef M. Gallmetzer).
- read whole integer tokens in selections ([488cda6](https://github.com/MolarVerse/PQAnalysis/commit/488cda6870ce4a82475b7d05274232b37ce58491) by Josef M. Gallmetzer).
- accept the documented numeric hessian_sign values ([7a35ebf](https://github.com/MolarVerse/PQAnalysis/commit/7a35ebf01798f40f380f58eaa85aa27078a6eaa9) by Josef M. Gallmetzer).
- preserve lazy analysis exports ([576eb98](https://github.com/MolarVerse/PQAnalysis/commit/576eb9835df09d292508f49f17dc38e7ea98abd7) by Josef M. Gallmetzer).
- handle empty files in MSD batches ([c5abf83](https://github.com/MolarVerse/PQAnalysis/commit/c5abf8329fe2b55617e4e0a0c7a9e648ed1c7f53) by Josef M. Gallmetzer).
- enforce momentum batch memory cap ([b2317d0](https://github.com/MolarVerse/PQAnalysis/commit/b2317d0ed5f53c01041cc1c8c0f6f7f81dffbf55) by Josef M. Gallmetzer).

### Performance Improvements

- reduce RDF and MSD trajectory setup ([b67b7fc](https://github.com/MolarVerse/PQAnalysis/commit/b67b7fc47e8b7b4cebd0ee4f4a15367bf514e7b0) by Josef M. Gallmetzer).
- reduce CLI startup overhead ([17a2bd2](https://github.com/MolarVerse/PQAnalysis/commit/17a2bd28dca2cfd2957a91ba6663ee116040f00d) by Josef M. Gallmetzer).
- batch exact RDF frames per worker ([2065bb8](https://github.com/MolarVerse/PQAnalysis/commit/2065bb82b45323f1330c72fa58ca826af60b5ef2) by Josef M. Gallmetzer).
- avoid repeated selection parser setup ([d2ac0c9](https://github.com/MolarVerse/PQAnalysis/commit/d2ac0c9a42afbde40b9de8b7f78e04b3125ef921) by Josef M. Gallmetzer).
- parallelize exact RDF frame histograms ([c81b5bb](https://github.com/MolarVerse/PQAnalysis/commit/c81b5bb7bce41bda16df9150cda4f5a361395772) by Josef M. Gallmetzer).
- batch exact MSD position loading ([645c927](https://github.com/MolarVerse/PQAnalysis/commit/645c9270b75fb7f9a2d098df17d353c8d778df1f) by Josef M. Gallmetzer).
- batch exact VACF velocity loading ([3022e26](https://github.com/MolarVerse/PQAnalysis/commit/3022e26ab6536e7cad8bf317efc17573f760a19c) by Josef M. Gallmetzer).
- fuse exact momentum parsing and reduction ([25063f2](https://github.com/MolarVerse/PQAnalysis/commit/25063f2f2016db8ab4f832da804f878525186986) by Josef M. Gallmetzer).
- add bounded trajectory batch parsing ([9e1d8b2](https://github.com/MolarVerse/PQAnalysis/commit/9e1d8b22052c67856b6733279c1f52f465667ac3) by Josef M. Gallmetzer).
- make momentum bitwise-exact with equipartition ([e8ee119](https://github.com/MolarVerse/PQAnalysis/commit/e8ee11941a98c1e5ed52996d1c890a8f1f9c51a2) by Josef M. Gallmetzer).
- make RDF bitwise-exact with legacy RDF ([8aa924c](https://github.com/MolarVerse/PQAnalysis/commit/8aa924cbea4bcb163a18a00c81940fa579bfad58) by Josef M. Gallmetzer).
- make MSD bitwise-exact with Diffcalc ([f78f0bd](https://github.com/MolarVerse/PQAnalysis/commit/f78f0bd7a39e3c4cc3e57012ffd01e7422e15a09) by Josef M. Gallmetzer).
- make VACF bitwise-exact with FreqCalc ([9c1addb](https://github.com/MolarVerse/PQAnalysis/commit/9c1addbf4cec88fe9f1422bc2bc6681b322e21d9) by Josef M. Gallmetzer).
- add direct float64 trajectory parsing ([bb4ea5f](https://github.com/MolarVerse/PQAnalysis/commit/bb4ea5fb3af96d7f390423010d9bc38eb7dc9cfb) by Josef M. Gallmetzer).

### Tests

- cover deferred CLI edge paths ([15aa5a2](https://github.com/MolarVerse/PQAnalysis/commit/15aa5a2e2800accd8fb7521e204cbdac4d2e57ad) by Josef M. Gallmetzer).
- cover batched analysis fallbacks ([4477efd](https://github.com/MolarVerse/PQAnalysis/commit/4477efda54c5ac8cadf2e70864661abce977fec8) by Josef M. Gallmetzer).

## [v1.5.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.5.0) - 2026-08-07

<small>[Compare with v1.4.1](https://github.com/MolarVerse/PQAnalysis/compare/v1.4.1...v1.5.0)</small>

### Features

- make XVG analysis exports reversible ([1b24af6](https://github.com/MolarVerse/PQAnalysis/commit/1b24af62b33511ff999af6ac02992d238b3ce123) by Josef M. Gallmetzer).
- add analysis table exports ([36c80ba](https://github.com/MolarVerse/PQAnalysis/commit/36c80bad0298557406b7b14b8f5e7e4f53391d52) by Josef M. Gallmetzer).

### Bug Fixes

- keep analysis tables readable ([55c9146](https://github.com/MolarVerse/PQAnalysis/commit/55c9146f62993425af950e11380897ea19defd40) by Josef M. Gallmetzer).
- preserve existing output filenames ([ed63dcd](https://github.com/MolarVerse/PQAnalysis/commit/ed63dcdf4520318fff2b12279e42187969ca2d9a) by Josef M. Gallmetzer).
- document analysis output columns ([92b9c25](https://github.com/MolarVerse/PQAnalysis/commit/92b9c25081579a77ac9be940bc1da53418f779f4) by Josef M. Gallmetzer).

### Style

- use Unicode analysis symbols ([27fc148](https://github.com/MolarVerse/PQAnalysis/commit/27fc148696b7cff200a202c48a5bbd3727b3ae7b) by Josef M. Gallmetzer).
- format analysis output headers ([3f98af4](https://github.com/MolarVerse/PQAnalysis/commit/3f98af4fab2f0b7c82e9b2af18e4d7a8c5554875) by Josef M. Gallmetzer).

## [v1.4.1](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.4.1) - 2026-07-30

<small>[Compare with v1.4.0](https://github.com/MolarVerse/PQAnalysis/compare/v1.4.0...v1.4.1)</small>

### Chore

- canonicalize contributor identity ([71dc6af](https://github.com/MolarVerse/PQAnalysis/commit/71dc6af9d9c20c106a877035d053a06f035acf41) by Josef M. Gallmetzer).

### Bug Fixes

- support single-entry PQ info rows (#160) ([ff24cc7](https://github.com/MolarVerse/PQAnalysis/commit/ff24cc7199e651ec913df13dacad434bb68f5fc2) by Josef M. Gallmetzer).

## [v1.4.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.4.0) - 2026-07-22

<small>[Compare with v1.3.0](https://github.com/MolarVerse/PQAnalysis/compare/v1.3.0...v1.4.0)</small>

### Build

- add conda-forge recipe (#152) ([2a7715e](https://github.com/MolarVerse/PQAnalysis/commit/2a7715e5accfe963c492863074fa53609e2d2e30) by Josef M. Gallmetzer).

### Continuous Integration

- update sigstore release action ([f8c15f6](https://github.com/MolarVerse/PQAnalysis/commit/f8c15f6f3e2d99e55c333da38951da3a3f4fa4f7) by Josef M. Gallmetzer).

### Docs

- fix sphinx build warnings (#149) ([ce0861c](https://github.com/MolarVerse/PQAnalysis/commit/ce0861c85cd288eac406716b81a27674c6d495d8) by Josef M. Gallmetzer).

### Features

- add PQ optimizer output reader (#156) ([bd653f9](https://github.com/MolarVerse/PQAnalysis/commit/bd653f9dbd2c8f5b0455ac8e664c20a441931b67) by Josef M. Gallmetzer).
- add MSD, VACF, spectrum and momentum analyses (faster than legacy C) (#153) ([f80bd9c](https://github.com/MolarVerse/PQAnalysis/commit/f80bd9c69595b69fbcf82fd6b091fa17732537e3) by Josef M. Gallmetzer).
- add vibrational analysis workflow (#150) ([58df167](https://github.com/MolarVerse/PQAnalysis/commit/58df167b57b897f8ad643cf7915dfa65e6ab78ba) by Josef M. Gallmetzer).

### Bug Fixes

- skip self pairs in RDF (#151) ([d5637f9](https://github.com/MolarVerse/PQAnalysis/commit/d5637f946b422839556602347416921a216a3217) by Josef M. Gallmetzer).
- support qmcfc input continuation (#148) ([a96859d](https://github.com/MolarVerse/PQAnalysis/commit/a96859decac7c42a290a3b67c6e5aee6bdbd0e38) by Josef M. Gallmetzer).

### Tests

- stabilize linear MSD fit tolerance (#158) ([a9d2922](https://github.com/MolarVerse/PQAnalysis/commit/a9d2922b599eedd49935f7dcd11af8fadef2d64d) by Josef M. Gallmetzer).

## [v1.3.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.3.0) - 2026-06-22

<small>[Compare with v1.2.5](https://github.com/MolarVerse/PQAnalysis/compare/v1.2.5...v1.3.0)</small>

### Continuous Integration

- refresh apt before docs build ([cf48554](https://github.com/MolarVerse/PQAnalysis/commit/cf48554da9212bfb1f9384b70fb5ad36ee6b9812) by Josef M. Gallmetzer).
- align dev with main ([51a1b1d](https://github.com/MolarVerse/PQAnalysis/commit/51a1b1d09151c5eea2e7277a81799e3a3b50e5a2) by Josef M. Gallmetzer).
- standardize pull request workflow (#129) ([dad5494](https://github.com/MolarVerse/PQAnalysis/commit/dad549423917daad5bdd7fe0e6c81f89e2bc1e75) by Josef M. Gallmetzer).
- standardize pull request workflow ([cf62ab4](https://github.com/MolarVerse/PQAnalysis/commit/cf62ab4d18e08471dc63648609f665f97b018cb9) by Josef M. Gallmetzer).

### Features

- add extxyz output profiles ([1c1bc0d](https://github.com/MolarVerse/PQAnalysis/commit/1c1bc0dbc49f3ae8b6013601fc3eb40fc1f17cd4) by Josef M. Gallmetzer).
- infer moltypes for xyz2rst (#141) ([709047f](https://github.com/MolarVerse/PQAnalysis/commit/709047fb7877edea4fddaa42d9f28dfe4e685ce5) by Josef M. Gallmetzer).
- convert pq trajectories to extxyz (#139) ([27f58ea](https://github.com/MolarVerse/PQAnalysis/commit/27f58eab2be40ba708e989f7a3375e1bc1b5bc84) by Josef M. Gallmetzer).
- read extended xyz trajectories (#138) ([bd7e4ae](https://github.com/MolarVerse/PQAnalysis/commit/bd7e4aeb8af2768576b4cfed2a3fe0bdcf9f9e45) by Josef M. Gallmetzer).
- continue inputs with unnumbered starts (#137) ([c1b7e14](https://github.com/MolarVerse/PQAnalysis/commit/c1b7e149326d0b8575bfa91cee175f642fbf83a1) by Josef M. Gallmetzer).
- add box file reader (#136) ([1bcfc69](https://github.com/MolarVerse/PQAnalysis/commit/1bcfc69792a78f3b4fcab2d56154484f27cb1201) by Josef M. Gallmetzer).
- add simulation time utility (#135) ([15b0b4b](https://github.com/MolarVerse/PQAnalysis/commit/15b0b4be42ab48ea5c40dc588c216bacfd0a4128) by Josef M. Gallmetzer).
- support multiple energy files (#132) ([f6d49c5](https://github.com/MolarVerse/PQAnalysis/commit/f6d49c5cdf125976879a6aa399031fb6af7e82a0) by Josef M. Gallmetzer).
- infer rdf topology files (#131) ([dfa4b8a](https://github.com/MolarVerse/PQAnalysis/commit/dfa4b8a7fd91d85f9a83873718638728f741562a) by Josef M. Gallmetzer).

### Bug Fixes

- preserve atom label case in gen conversions (#133) ([0a95159](https://github.com/MolarVerse/PQAnalysis/commit/0a951591732a6c337a536ecb5512fea353cf3344) by Josef M. Gallmetzer).
- resolve rdf input handling (#114) ([4898f8f](https://github.com/MolarVerse/PQAnalysis/commit/4898f8fef9a11fd2ad62c98510354f0c06c17161) by Stefanie Kröll).

### Code Refactoring

- split trajectory frame reader (#127) ([ff5d156](https://github.com/MolarVerse/PQAnalysis/commit/ff5d156b0f9002a9e1221765a8aa9bb7354960a6) by Josef M. Gallmetzer).

### Tests

- cover release type-checking output (#140) ([272d315](https://github.com/MolarVerse/PQAnalysis/commit/272d315a7fcf8b51ad1e159c68fb136fee9d2e48) by Josef M. Gallmetzer).

## [v1.2.5](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.2.5) - 2025-02-12

<small>[Compare with v1.2.4](https://github.com/MolarVerse/PQAnalysis/compare/v1.2.4...v1.2.5)</small>

## [v1.2.4](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.2.4) - 2025-02-12

<small>[Compare with v1.2.3](https://github.com/MolarVerse/PQAnalysis/compare/v1.2.3...v1.2.4)</small>

## [v1.2.3](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.2.3) - 2025-02-12

<small>[Compare with v1.2.2](https://github.com/MolarVerse/PQAnalysis/compare/v1.2.2...v1.2.3)</small>

### Docs

- deleted no-box from output docs ([75eaea6](https://github.com/MolarVerse/PQAnalysis/commit/75eaea696a6c7749ed4c77613c5486afa64a5d5d) by Josef M. Gallmetzer).

### Features

- add randomization option for atom positions in xyz2rst conversion ([1a1ddc5](https://github.com/MolarVerse/PQAnalysis/commit/1a1ddc5a4e5cd92e54a7d31f37175beced726d57) by Josef M. Gallmetzer).
- add xyz2rst conversion functionality and related tests ([598473f](https://github.com/MolarVerse/PQAnalysis/commit/598473f778f4d8bc5cbbd22de43d64f2ea1d7a23) by Josef M. Gallmetzer).
- add API function to write restart files with specified format and mode ([184061e](https://github.com/MolarVerse/PQAnalysis/commit/184061ed0d51e01bd6de8d04bf8c441e900e6cf2) by Josef M. Gallmetzer).

### Bug Fixes

- update restart writer to conditionally include velocities and forces in output ([bd841d9](https://github.com/MolarVerse/PQAnalysis/commit/bd841d9999b44a6b183e0ccabc5be718038ef3f5) by Josef M. Gallmetzer).
- prevent adding box line for vacuum cells in restart file writer ([5f72088](https://github.com/MolarVerse/PQAnalysis/commit/5f72088c4179d6d68d7b62fa2fc4eccc8d31c0af) by Josef M. Gallmetzer).

### Tests

- add test for XYZ2RstCLI program name ([e7ee687](https://github.com/MolarVerse/PQAnalysis/commit/e7ee687425532e264d47526939c19839c4d9e868) by Josef M. Gallmetzer).
- add test for xyz2rst conversion with no box data ([119cd75](https://github.com/MolarVerse/PQAnalysis/commit/119cd75829abe883b1bde55982d03975337bade4) by Josef M. Gallmetzer).

## [v1.2.2](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.2.2) - 2025-01-07

<small>[Compare with v1.2.1](https://github.com/MolarVerse/PQAnalysis/compare/v1.2.1...v1.2.2)</small>

### Bug Fixes

- pin setuptools to version 70.0.0 so that pytest.sh works ([1b191f8](https://github.com/MolarVerse/PQAnalysis/commit/1b191f85f70d6f1b1eff0e79b600fba454527cc4) by Josef M. Gallmetzer).
- update setuptools dependency to latest version ([d52e3b1](https://github.com/MolarVerse/PQAnalysis/commit/d52e3b16475a5cb8a48dab4fd51f79987d0b25f7) by Josef M. Gallmetzer).
- Replace np.in1d with np.isin for improved index matching ([cfb5391](https://github.com/MolarVerse/PQAnalysis/commit/cfb5391f76fd7fee0d31b32c8f0b6d4d41455036) by Josef M. Gallmetzer).
- Disable iteration over AtomicSystem object ([bfb0365](https://github.com/MolarVerse/PQAnalysis/commit/bfb03650c9d902569b20409bae534b4a3fc3189f) by Josef M. Gallmetzer).

## [v1.2.1](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.2.1) - 2024-08-19

<small>[Compare with v1.1.2](https://github.com/MolarVerse/PQAnalysis/compare/v1.1.2...v1.2.1)</small>

### Continuous Integration

- changelog creation for release notes refactored - hopefully working now ([66c8d81](https://github.com/MolarVerse/PQAnalysis/commit/66c8d81326284f56bb691e4ca93715a2a0feb925) by 97gamjak).
- release.yml refactored to not upload wheels due to cython build ([0f42fe5](https://github.com/MolarVerse/PQAnalysis/commit/0f42fe56a45dc9523b6fa10df7052d036e1f6cdb) by 97gamjak).
- docs deployment should now work - cleaned up ([b7b318c](https://github.com/MolarVerse/PQAnalysis/commit/b7b318cac046df71d3311ea8fa9277bcea89cf71) by 97gamjak).
- release deployment should work now with cython ([824fbfe](https://github.com/MolarVerse/PQAnalysis/commit/824fbfe68348496c31b34236cb40b2050a3c1ba1) by 97gamjak).
- docs deployment fix to wotk now with cython ([a91d409](https://github.com/MolarVerse/PQAnalysis/commit/a91d40920d7e1d7306f8827eaa24e18c57c6edee) by 97gamjak).
- added optional dependency of setuptools to [test] ([93f92de](https://github.com/MolarVerse/PQAnalysis/commit/93f92dedc45f1966ba81a889541160b13365834f) by Jakob Gamper).

### Bug Fixes

- small bugfix in release.yml ([2283160](https://github.com/MolarVerse/PQAnalysis/commit/22831600199e05e9ebbad3df339b9e52e37d216c) by Jakob Gamper).

### Performance Improvements

- Improved performance of reading a trajectory ([d929e0b](https://github.com/MolarVerse/PQAnalysis/commit/d929e0b4616632c7a1a00dd0bba659262f1a3be1) by Jakob Gamper).

### Tests

- fixed pytest.sh ([64967cc](https://github.com/MolarVerse/PQAnalysis/commit/64967ccb5bd78929ea3b27835bf7d8fe8d4bf771) by Jakob Gamper).
- pytest.sh updated to work with cython files ([9950f34](https://github.com/MolarVerse/PQAnalysis/commit/9950f34c9164d47349a5d5dd158bbc04506f2e45) by Jakob Gamper).

## [v1.1.2](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.1.2) - 2024-06-07

<small>[Compare with v1.1.1](https://github.com/MolarVerse/PQAnalysis/compare/v1.1.1...v1.1.2)</small>

### Continuous Integration

- updated release.yml ([6d51565](https://github.com/MolarVerse/PQAnalysis/commit/6d515650ddc42a7ad6e91bd364c4540dc39697b3) by Jakob Gamper).

## [v1.1.1](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.1.1) - 2024-06-05

<small>[Compare with v1.0.12](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.12...v1.1.1)</small>

### Docs

- Add logo to README.md ([a9aef5c](https://github.com/MolarVerse/PQAnalysis/commit/a9aef5c7afaf9a385c12f477d8b3aba66c41b7bb) by Josef M. Gallmetzer).

### Features

- Added possibility to import directly from PQAnalysis.topology.bonded_topology ([75de1a7](https://github.com/MolarVerse/PQAnalysis/commit/75de1a79eec3837cb3a89871d1ee48acf0399da3) by Jakob Gamper).
- merging topologies keeps now comments of data line to be more consistent with input topologies ([d224e78](https://github.com/MolarVerse/PQAnalysis/commit/d224e789f2e82993a136c585c84cf019670bfe44) by Jakob Gamper).
- added utils function to check if line is a comment_line ([9852e33](https://github.com/MolarVerse/PQAnalysis/commit/9852e33dcda8042d0e3aef7ec726ef2f859ee6fe) by Jakob Gamper).

### Bug Fixes

- read_trajectory with constant topology did actually not use a constant topology approach ([ef1c5db](https://github.com/MolarVerse/PQAnalysis/commit/ef1c5dba3752f9308685e2ce57cb017e54e719e1) by Jakob Gamper).
- added missing linker output for linker bonds, angles, ... ([ad47473](https://github.com/MolarVerse/PQAnalysis/commit/ad47473f1e0f25374fc526eeb9e747108a9dc6a0) by Jakob Gamper).
- included possibility to write topology to stdout ([3ea89c3](https://github.com/MolarVerse/PQAnalysis/commit/3ea89c3213ef88b8d9462166e9de60147772928b) by Jakob Gamper).

### Tests

- added some unit tests for TopologyFileWriter ([65891fa](https://github.com/MolarVerse/PQAnalysis/commit/65891fa2138861024fbb2589cde061313187e7ee) by Jakob Gamper).
- first integration tests of add_molecules added ([bc9e694](https://github.com/MolarVerse/PQAnalysis/commit/bc9e694dc6e48dd31fb9cf1973354356904c53ae) by Jakob Gamper).
- added constant seed strategy for executing tests ([b6a14b7](https://github.com/MolarVerse/PQAnalysis/commit/b6a14b72e5a014b17047268da7a117e13e413a2a) by Jakob Gamper).
- added basic structure for integration tests ([729334e](https://github.com/MolarVerse/PQAnalysis/commit/729334e3b647f40937e17572effa79aec1a726c8) by Jakob Gamper).


## [v1.0.12](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.12) - 2024-06-01

<small>[Compare with v1.0.11](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.11...v1.0.12)</small>

### Continuous Integration

- cleaned up release.yml for changelog writing ([f9435b4](https://github.com/MolarVerse/PQAnalysis/commit/f9435b4dcc316e12e5f0a40a21510bb72fc1cc8d) by Jakob Gamper).

## [v1.0.11](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.11) - 2024-06-01

<small>[Compare with v1.0.10](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.10...v1.0.11)</small>

## [v1.0.10](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.10) - 2024-05-31

<small>[Compare with v1.0.9](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.9...v1.0.10)</small>

### Bug Fixes

- logger now working again also for cli tools ([e101a7b](https://github.com/MolarVerse/PQAnalysis/commit/e101a7b7d44e1ac41de1ee33c3058a22960ad33b) by Jakob Gamper).

## [v1.0.9](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.9) - 2024-05-31

<small>[Compare with v1.0.7](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.7...v1.0.9)</small>

### Code Refactoring

- Update release workflow to include PyPI publishing on tag pushes ([5233a6c](https://github.com/MolarVerse/PQAnalysis/commit/5233a6c6056099efa275493fc38cf986981d511c) by Jakob Gamper).
- Update release workflow to include CHANGELOG.md generation and commit ([b52c958](https://github.com/MolarVerse/PQAnalysis/commit/b52c958a10c8346022173d2ad989a7401d16f8aa) by Jakob Gamper).
- Update release workflow to include permissions and branch filtering ([ccfe903](https://github.com/MolarVerse/PQAnalysis/commit/ccfe903571a50f5cc6a2d1a693185c00c5b03810) by Jakob Gamper).
- Update release workflow to include CHANGELOG.md generation and commit [skip ci] ([1e81a5f](https://github.com/MolarVerse/PQAnalysis/commit/1e81a5f751e633a2d2ee50075104548580098839) by Jakob Gamper).

## [v1.0.7](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.7) - 2024-05-31

<small>[Compare with v1.0.6](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.6...v1.0.7)</small>

## [v1.0.6](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.6) - 2024-05-30

<small>[Compare with v1.0.5](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.5...v1.0.6)</small>

### Features

- Use vectorized allclose function for element-wise comparison of arrays in Cell class ([6e98813](https://github.com/MolarVerse/PQAnalysis/commit/6e9881395f00749d06e8e295b1f6d62c4161908d) by Josef M. Gallmetzer).
- Update isclose method to accept any object for comparison ([28b3bee](https://github.com/MolarVerse/PQAnalysis/commit/28b3bee6b202a7da2f3c90b0f8ec305c3fbbf8c7) by Jakob Gamper).
- Add pytest marker for utils module tests ([1100b4d](https://github.com/MolarVerse/PQAnalysis/commit/1100b4d00442285d4fb76eff5bd3ffc81792bbf2) by Jakob Gamper).
- Add vectorized allclose function for element-wise comparison of numpy arrays ([b70f7ef](https://github.com/MolarVerse/PQAnalysis/commit/b70f7efe0888140a1b8524d848eda0bb77d4ed2a) by Jakob Gamper).
- Add isclose method to Trajectory class for comparing trajectories ([22b5ff2](https://github.com/MolarVerse/PQAnalysis/commit/22b5ff2b8381f993291226a048bc632f98934069) by Josef M. Galletzer).
- Add functionality to add line comments to the topology ([041a5fc](https://github.com/MolarVerse/PQAnalysis/commit/041a5fcb0674c8aaa35d4d459cf5f080426b916e) by Jakob Gamper).

### Code Refactoring

- Update .github/workflows/pylint.yml to include .github/.pylint_cache in the commit ([15a5cfc](https://github.com/MolarVerse/PQAnalysis/commit/15a5cfcb7c65d717b4544e185158af76ae0220a9) by Jakob Gamper).
- Update isclose method to use smaller default tolerances for element-wise comparison ([48668aa](https://github.com/MolarVerse/PQAnalysis/commit/48668aa3f9783e675c248d4e9b8bfbd24f68469e) by Josef M. Gallmetzer).
- Update Trajectory class to use vectorized allclose function for element-wise comparison of arrays ([f9fa8c8](https://github.com/MolarVerse/PQAnalysis/commit/f9fa8c83e2c0a876d632aa59d8f4eb57b6a9587d) by Josef M. Gallmetzer).
- Use vectorized allclose function for element-wise comparison of arrays ([84608c4](https://github.com/MolarVerse/PQAnalysis/commit/84608c4421004c9476b559fd6265d3498ff7e3e9) by Josef M. Gallmetzer).

## [v1.0.5](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.5) - 2024-05-26

<small>[Compare with v1.0.4](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.4...v1.0.5)</small>

## [v1.0.4](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.4) - 2024-05-25

<small>[Compare with v1.0.3](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.3...v1.0.4)</small>

## [v1.0.3](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.3) - 2024-05-23

<small>[Compare with v1.0.2](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.2...v1.0.3)</small>

### Features

- Add benchmark for reading trajectories with different frame counts ([ec82fd6](https://github.com/MolarVerse/PQAnalysis/commit/ec82fd680dc0b6da8420ffff0df0834a935aad3a) by Josef M. Galletzer).

## [v1.0.2](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.2) - 2024-05-13

<small>[Compare with v1.0.1](https://github.com/MolarVerse/PQAnalysis/compare/v1.0.1...v1.0.2)</small>

## [v1.0.1](https://github.com/MolarVerse/PQAnalysis/releases/tag/v1.0.1) - 2024-05-13

<small>[Compare with v0.5.2](https://github.com/MolarVerse/PQAnalysis/compare/v0.5.2...v1.0.1)</small>

### Features

- Refactor TrajectoryReader window generator and test final ([7c36182](https://github.com/MolarVerse/PQAnalysis/commit/7c3618258e55103af2127d82bef31d1808f0619a) by Josef M. Gallmetzer).
- Improve error handling in TrajectoryReader class with tests (not complete) ([4f93e9a](https://github.com/MolarVerse/PQAnalysis/commit/4f93e9a72d967eeb50aaad6b6ef8ff359107b832) by Josef M. Gallmetzer).
- Add pop method to Trajectory class ([1a82dc3](https://github.com/MolarVerse/PQAnalysis/commit/1a82dc3b1c7d199d5ae1d08a10694f2ce654e1ee) by Josef M. Gallmetzer).
- Update Trajectory Reader class to include a window generator ([6de6d9d](https://github.com/MolarVerse/PQAnalysis/commit/6de6d9dba052f0fa960ee88af74481f6436fa560) by Josef M. Gallmetzer).
- add window_generator method to TrajectoryReader class ([95d0191](https://github.com/MolarVerse/PQAnalysis/commit/95d01910a038b802f3d262320508153556808ef3) by Josef M. Gallmetzer).
- Add custom_exception attribute to log records in CustomLogger class ([a0c9ce7](https://github.com/MolarVerse/PQAnalysis/commit/a0c9ce72ce524819734d01fbbc28dba14bac87c8) by Josef M. Gallmetzer).
- Update Trajectory class window method documentation ([5206384](https://github.com/MolarVerse/PQAnalysis/commit/5206384d111ca02e1a9cd568f5d03513036b40b4) by Josef M. Gallmetzer).
- Refactor window method in Trajectory class and adjusted the tests ([628e6bb](https://github.com/MolarVerse/PQAnalysis/commit/628e6bb808094ed9a2dffe1e7e8325a26eebba17) by Josef M. Galletzer).
- add window method to Trajectory class ([84bc46c](https://github.com/MolarVerse/PQAnalysis/commit/84bc46c34d981b22b1928e5a84867556e5293e11) by Josef M. Galletzer).

### Code Refactoring

- Initialize length_of_traj to 0 in TrajectoryReader constructor and tested file change ([f5834e1](https://github.com/MolarVerse/PQAnalysis/commit/f5834e19e13294712fc5747e4fbbefe798cb5831) by Josef M. Galletzer).

## [v0.5.2](https://github.com/MolarVerse/PQAnalysis/releases/tag/v0.5.2) - 2023-12-09

<small>[Compare with v0.5.1](https://github.com/MolarVerse/PQAnalysis/compare/v0.5.1...v0.5.2)</small>

## [v0.5.1](https://github.com/MolarVerse/PQAnalysis/releases/tag/v0.5.1) - 2023-11-28

<small>[Compare with v0.5.0](https://github.com/MolarVerse/PQAnalysis/compare/v0.5.0...v0.5.1)</small>

## [v0.5.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v0.5.0) - 2023-11-28

<small>[Compare with v0.4.0](https://github.com/MolarVerse/PQAnalysis/compare/v0.4.0...v0.5.0)</small>

## [v0.4.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v0.4.0) - 2023-11-12

<small>[Compare with v0.3.2](https://github.com/MolarVerse/PQAnalysis/compare/v0.3.2...v0.4.0)</small>

## [v0.3.2](https://github.com/MolarVerse/PQAnalysis/releases/tag/v0.3.2) - 2023-11-10

<small>[Compare with v0.3.1](https://github.com/MolarVerse/PQAnalysis/compare/v0.3.1...v0.3.2)</small>

## [v0.3.1](https://github.com/MolarVerse/PQAnalysis/releases/tag/v0.3.1) - 2023-11-09

<small>[Compare with v0.3.0](https://github.com/MolarVerse/PQAnalysis/compare/v0.3.0...v0.3.1)</small>

## [v0.3.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v0.3.0) - 2023-11-09

<small>[Compare with v0.2.0](https://github.com/MolarVerse/PQAnalysis/compare/v0.2.0...v0.3.0)</small>

## [v0.2.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v0.2.0) - 2023-10-31

<small>[Compare with v0.1.0](https://github.com/MolarVerse/PQAnalysis/compare/v0.1.0...v0.2.0)</small>

## [v0.1.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v0.1.0) - 2023-10-24

<small>[Compare with v0.0.0](https://github.com/MolarVerse/PQAnalysis/compare/v0.0.0...v0.1.0)</small>

## [v0.0.0](https://github.com/MolarVerse/PQAnalysis/releases/tag/v0.0.0) - 2023-10-23

<small>[Compare with first commit](https://github.com/MolarVerse/PQAnalysis/compare/e5b4d04ce4e5a3c6e910f027a1f443cb0fc1fb39...v0.0.0)</small>
