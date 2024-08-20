# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

<!-- insertion marker -->
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
