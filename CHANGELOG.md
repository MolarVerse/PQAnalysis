# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.0.6.1] - 2024-05-31
### :sparkles: New Features
- [`84bc46c`](https://github.com/MolarVerse/PQAnalysis/commit/84bc46c34d981b22b1928e5a84867556e5293e11) - add window method to Trajectory class *(commit by [@galjos](https://github.com/galjos))*
- [`628e6bb`](https://github.com/MolarVerse/PQAnalysis/commit/628e6bb808094ed9a2dffe1e7e8325a26eebba17) - Refactor window method in Trajectory class and adjusted the tests *(commit by [@galjos](https://github.com/galjos))*
- [`5206384`](https://github.com/MolarVerse/PQAnalysis/commit/5206384d111ca02e1a9cd568f5d03513036b40b4) - Update Trajectory class window method documentation *(commit by [@galjos](https://github.com/galjos))*
- [`a0c9ce7`](https://github.com/MolarVerse/PQAnalysis/commit/a0c9ce72ce524819734d01fbbc28dba14bac87c8) - Add custom_exception attribute to log records in CustomLogger class *(commit by [@galjos](https://github.com/galjos))*
- [`95d0191`](https://github.com/MolarVerse/PQAnalysis/commit/95d01910a038b802f3d262320508153556808ef3) - add window_generator method to TrajectoryReader class *(commit by [@galjos](https://github.com/galjos))*
- [`6de6d9d`](https://github.com/MolarVerse/PQAnalysis/commit/6de6d9dba052f0fa960ee88af74481f6436fa560) - Update Trajectory Reader class to include a window generator *(commit by [@galjos](https://github.com/galjos))*
- [`1a82dc3`](https://github.com/MolarVerse/PQAnalysis/commit/1a82dc3b1c7d199d5ae1d08a10694f2ce654e1ee) - Add pop method to Trajectory class *(commit by [@galjos](https://github.com/galjos))*
- [`4f93e9a`](https://github.com/MolarVerse/PQAnalysis/commit/4f93e9a72d967eeb50aaad6b6ef8ff359107b832) - Improve error handling in TrajectoryReader class with tests (not complete) *(commit by [@galjos](https://github.com/galjos))*
- [`7c36182`](https://github.com/MolarVerse/PQAnalysis/commit/7c3618258e55103af2127d82bef31d1808f0619a) - Refactor TrajectoryReader window generator and test final *(commit by [@galjos](https://github.com/galjos))*
- [`ec82fd6`](https://github.com/MolarVerse/PQAnalysis/commit/ec82fd680dc0b6da8420ffff0df0834a935aad3a) - Add benchmark for reading trajectories with different frame counts *(commit by [@galjos](https://github.com/galjos))*
- [`041a5fc`](https://github.com/MolarVerse/PQAnalysis/commit/041a5fcb0674c8aaa35d4d459cf5f080426b916e) - Add functionality to add line comments to the topology *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`22b5ff2`](https://github.com/MolarVerse/PQAnalysis/commit/22b5ff2b8381f993291226a048bc632f98934069) - Add isclose method to Trajectory class for comparing trajectories *(commit by [@galjos](https://github.com/galjos))*
- [`b70f7ef`](https://github.com/MolarVerse/PQAnalysis/commit/b70f7efe0888140a1b8524d848eda0bb77d4ed2a) - Add vectorized allclose function for element-wise comparison of numpy arrays *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`1100b4d`](https://github.com/MolarVerse/PQAnalysis/commit/1100b4d00442285d4fb76eff5bd3ffc81792bbf2) - Add pytest marker for utils module tests *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`28b3bee`](https://github.com/MolarVerse/PQAnalysis/commit/28b3bee6b202a7da2f3c90b0f8ec305c3fbbf8c7) - Update isclose method to accept any object for comparison *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`6e98813`](https://github.com/MolarVerse/PQAnalysis/commit/6e9881395f00749d06e8e295b1f6d62c4161908d) - Use vectorized allclose function for element-wise comparison of arrays in Cell class *(commit by [@galjos](https://github.com/galjos))*

### :bug: Bug Fixes
- [`0594013`](https://github.com/MolarVerse/PQAnalysis/commit/0594013961e2db7404f71f4bcee2e81b14a1fa75) - energy setattr class to object *(commit by [@galjos](https://github.com/galjos))*

### :recycle: Refactors
- [`f5834e1`](https://github.com/MolarVerse/PQAnalysis/commit/f5834e19e13294712fc5747e4fbbefe798cb5831) - Initialize length_of_traj to 0 in TrajectoryReader constructor and tested file change *(commit by [@galjos](https://github.com/galjos))*
- [`84608c4`](https://github.com/MolarVerse/PQAnalysis/commit/84608c4421004c9476b559fd6265d3498ff7e3e9) - Use vectorized allclose function for element-wise comparison of arrays *(commit by [@galjos](https://github.com/galjos))*
- [`f9fa8c8`](https://github.com/MolarVerse/PQAnalysis/commit/f9fa8c83e2c0a876d632aa59d8f4eb57b6a9587d) - Update Trajectory class to use vectorized allclose function for element-wise comparison of arrays *(commit by [@galjos](https://github.com/galjos))*
- [`43ab020`](https://github.com/MolarVerse/PQAnalysis/commit/43ab020c2e379daa7421645986798f24b5e9f9b0) - Update isclose method to use smaller default tolerances for element-wise comparison *(commit by [@galjos](https://github.com/galjos))*
- [`48668aa`](https://github.com/MolarVerse/PQAnalysis/commit/48668aa3f9783e675c248d4e9b8bfbd24f68469e) - Update isclose method to use smaller default tolerances for element-wise comparison *(commit by [@galjos](https://github.com/galjos))*
- [`15a5cfc`](https://github.com/MolarVerse/PQAnalysis/commit/15a5cfcb7c65d717b4544e185158af76ae0220a9) - Update .github/workflows/pylint.yml to include .github/.pylint_cache in the commit *(commit by [@97gamjak](https://github.com/97gamjak))*

### :white_check_mark: Tests
- [`97a29f5`](https://github.com/MolarVerse/PQAnalysis/commit/97a29f55825531c41391bd995ad00c7486211b30) - added test cases for cell reader *(commit by [@galjos](https://github.com/galjos))*

### :wrench: Chores
- [`37bfdf6`](https://github.com/MolarVerse/PQAnalysis/commit/37bfdf60f02b5e2118b83c756c0f255bd226b5cb) - Improve error handling in Trajectory class window method with logger *(commit by [@galjos](https://github.com/galjos))*
- [`f236544`](https://github.com/MolarVerse/PQAnalysis/commit/f236544059f5dcaa1e70405f836bae64dbc15be0) - Update Trajectory class window parameters names in trajectory start and trajectory stop *(commit by [@galjos](https://github.com/galjos))*
- [`276dc8c`](https://github.com/MolarVerse/PQAnalysis/commit/276dc8cf585a31e84603cf3735ea79d0f011e02f) - Add 'six' package to dependencies *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`c39cdf0`](https://github.com/MolarVerse/PQAnalysis/commit/c39cdf03432a9fb6f683176711ca98802700d12a) - Refactor trajectory_reader.py for improved readability and maintainability *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`955c360`](https://github.com/MolarVerse/PQAnalysis/commit/955c360ea2e78af997f344fd0071fafd6fee92f8) - Update codecov-action to v4 in CI workflow *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`f7bb2e8`](https://github.com/MolarVerse/PQAnalysis/commit/f7bb2e8b1884d55b754427f6bdc90c86c6ccf211) - Update permissions in CI workflow *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`a461a19`](https://github.com/MolarVerse/PQAnalysis/commit/a461a19bdc76af807de3972313856413caf6db58) - Update permissions in CI workflow *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`f556c95`](https://github.com/MolarVerse/PQAnalysis/commit/f556c9581d0781e02df43d984e3064a3b15bfdb0) - Add `isclose` method to Cell class for comparing cell objects *(commit by [@galjos](https://github.com/galjos))*
- [`26ede9d`](https://github.com/MolarVerse/PQAnalysis/commit/26ede9d3301b114c8888cda04cc2eeb8a99d6efb) - Add isclose method to AtomicSystem class for comparing AtomicSystem objects *(commit by [@galjos](https://github.com/galjos))*
- [`088c7bc`](https://github.com/MolarVerse/PQAnalysis/commit/088c7bc9ed1ee7a24170f4b78414dde30f7503ec) - Add unit test for allclose_vectorized function in test_math.py *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`f21e07d`](https://github.com/MolarVerse/PQAnalysis/commit/f21e07d5f4050c9093aa011c3742ed35bed4b114) - Add utils marker to pytest.ini *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`4dbe4a0`](https://github.com/MolarVerse/PQAnalysis/commit/4dbe4a0267381b5835957efde141d89530110795) - Add Bool as new type for beartype and Format type hints in types.py for better readability *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`aeb9acf`](https://github.com/MolarVerse/PQAnalysis/commit/aeb9acfadcd018445111a22ed6471264a61d7b3a) - Add Bool type to AtomicSystem and Trajectory classes for improved comparison functionality *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`834bedc`](https://github.com/MolarVerse/PQAnalysis/commit/834bedc3c445d90e80271188e61d53bdff89f924) - bugfix in Bool type definition *(commit by [@97gamjak](https://github.com/97gamjak))*
- [`7b4829f`](https://github.com/MolarVerse/PQAnalysis/commit/7b4829f1d9c3ab6a7ca2b6506455040dc089a179) - Add error handling for pytest.sh script *(commit by [@97gamjak](https://github.com/97gamjak))*

[v1.0.6.1]: https://github.com/MolarVerse/PQAnalysis/compare/v0.0.0...v1.0.6.1
