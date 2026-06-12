---
name: Release Process
description: Guide for performing the DirectXMath release process. Use this skill when asked to help with releasing a new version, publishing packages, or updating ports.
---

# Release Process

## Prerequisites

- All changes merged into the `main` branch with all tests passing.
- GPG signing configured for your GitHub account (for verified tags).
- Access to the MSCodeHub mirror repository and Azure DevOps pipelines.
- Local repository:
  - VCPKG at `d:\vcpkg` (synced with `main` branch)
- PATs will be needed for scripts that access GitHub and ADO.

<!-- markdownlint-disable MD029 -->
## Steps

### Phase 1: Prepare Release

1. Git pull the local repository to ensure it is up to date with the `main` branch.
2. Run the PowerShell script `build\preparerelease.ps1` which will generate a topic branch for the release, update the version number in `CMakeLists.txt`, the `README.md` file, the release notes in the nuspec files, and create a stub in the `CHANGELOG.md` file for the new release.
3. Edit the `CHANGELOG.md` file to update it with a summary of changes.
4. Submit the topic branch for review and merge into `main` once approved. Allow the GitHub Actions workflows and the Azure DevOps pipelines to complete successfully before proceeding.

### Phase 2: Tag and Create GitHub Release

5. Run the PowerShell script `build\completerelease.ps1` which will set a tag on the project repo and the test repo, and create a release on GitHub with the release notes from `CHANGELOG.md`. Ensure you have set up GPG signing for your GitHub account so that the tags will be verified.
6. Git pull the local repository to ensure it is up to date with the `main` branch. Be sure to include `--tags`.

### Phase 3: MSCodeHub and Signed Binaries

7. Push the `main` branch to the MSCodeHub mirror repository. Be sure to include `--tags`.
8. Create a PR on MSCodeHub from the `main` branch to the `release` branch.
9. Merge the PR on MSCodeHub to update the release branch, which will trigger the Azure DevOps pipeline to build the NuGet package.

### Phase 4: Source Archive Signing

10. Download the GitHub source .zip archive from the release. Unzip and compare to the local repo to ensure it matches — keep in mind there may be some CR/LF differences.
11. Run minisign on the .zip to generate a signature file, and upload the signature file to the release assets.

### Phase 5: NuGet Validation and Publishing

12. Download the NuGet package from the ADO pipeline build, and review it with a local project or using NuGet Package Explorer.
13. Upload the nupkg to nuget.org via the website.

### Phase 6: VCPKG Port Update

14. Git pull a local repository of VCPKG to `d:\vcpkg` in sync with the `main` branch of the VCPKG repository.
15. Run the PowerShell script `build\updatevcpkg.ps1` to update the DirectXMath port in VCPKG with the new release version. This will edit the files in `ports\directxmath`.
16. Test the VCPKG port using the script at `assets/vcpkgdxmath.cmd` (in this skill folder). Copy it to `d:\vcpkg` and run from there after bootstrapping VCPKG.
17. Run `.\vcpkg x-add-version directxmath` to update the VCPKG versioning history.
18. Submit a PR to the VCPKG GitHub repository to update the DirectXMath port. The PR will be reviewed and merged by the VCPKG maintainers.

### Phase 7: Windows SDK Update

19. For the DirectXMath release to be included in the next Windows SDK, prepare a PR for the MSCodeHub project from the `main` branch to the `ms_sdk_release` branch. When the PR is complete, the Azure DevOps pipeline will automatically build vpack and submit a PR for further review.

### Phase 8: Finalize

20. When fully completed, be sure to update the GitHub release with links to the matching NuGet packages and the VCPKG port.
21. Contact the Visual C++ compiler team and have them update `src/qa/suites/common/TestAssets.xml` to use the release's commit ids for DirectXMath and the Test Suite to be included in their test coverage.

## Key Scripts

| Script | Purpose |
| --- | --- |
| `build\preparerelease.ps1` | Creates topic branch, updates version numbers and changelog stub |
| `build\completerelease.ps1` | Sets tags, creates GitHub release from changelog |
| `build\updatevcpkg.ps1` | Updates DirectXMath VCPKG port files |
| `assets\vcpkgdxmath.cmd` | Tests VCPKG port across all triplets and features |
