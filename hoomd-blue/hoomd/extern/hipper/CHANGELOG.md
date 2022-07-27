# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2020-10-16
### Added
- Wrapper around CUB library, which can be accessed like `hipper::cub`.
### Changed
- C++14 standard is now required.
### Fixed
- Inline template specifications for `threadRang`.

## [0.1.2] - 2020-07-01
### Fixed
- Remove `BUILD_TESTING` from CMake setup.
- Add missing changelog entries for [0.1.1].

## [0.1.1] - 2020-06-30
### Added
- Testing framework for CUDA and HIP targets.
### Fixed
- Definition of `limit_t` in HIP.
- Example in README.

## [0.1.0] - 2020-05-19
### Added
- Basic runtime wrappers for CUDA and HIP.

[Unreleased]: https://github.com/mphowardlab/hipper/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/mphowardlab/hipper/releases/tag/v0.2.0
[0.1.2]: https://github.com/mphowardlab/hipper/releases/tag/v0.1.2
[0.1.1]: https://github.com/mphowardlab/hipper/releases/tag/v0.1.1
[0.1.0]: https://github.com/mphowardlab/hipper/releases/tag/v0.1.0
