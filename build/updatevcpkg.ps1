<#

.NOTES
Copyright (c) Microsoft Corporation.
Licensed under the MIT License.

.SYNOPSIS
Updates the vcpkg port for DirectXMath to match a GitHub release.

.DESCRIPTION
This script updates the vcpkg port at D:\vcpkg\ports\directxmath to match a specific
GitHub release by tag. It updates the version-date in vcpkg.json and the tag and SHA512
hash in portfile.cmake for the source archive.

.PARAMETER Tag
The GitHub release tag (e.g., 'may2026', 'mar2026'). Defaults to the latest tag.

.LINK
https://github.com/microsoft/DirectXMath/wiki

#>

param(
    [string]$Tag = ""
)

$repoRoot = Split-Path -Path $PSScriptRoot -Parent
$portDir = "D:\vcpkg\ports\directxmath"
$vcpkgJson = Join-Path $portDir "vcpkg.json"
$portfile = Join-Path $portDir "portfile.cmake"

if ((-Not (Test-Path $vcpkgJson)) -Or (-Not (Test-Path $portfile))) {
    Write-Error "ERROR: Cannot find vcpkg port files at $portDir" -ErrorAction Stop
}

# Determine tag from latest git tag if not provided
if ($Tag.Length -eq 0) {
    $Tag = (git --no-pager -C $repoRoot tag --sort=-creatordate | Select-Object -First 1).Trim()
    if ($Tag.Length -eq 0) {
        Write-Error "ERROR: Failed to determine latest tag!" -ErrorAction Stop
    }
}

Write-Host "Release Tag: $Tag"

# Get version date from the git tag date
$tagDateStr = (git --no-pager -C $repoRoot log -1 --format=%ai $Tag).Trim()
if ([string]::IsNullOrEmpty($tagDateStr)) {
    Write-Error "ERROR: Failed to get date for tag $Tag!" -ErrorAction Stop
}

$versionDate = ([datetime]::Parse($tagDateStr)).ToString("yyyy-MM-dd")

Write-Host "Version Date: $versionDate"

# --- Update vcpkg.json ---
Write-Host "`nUpdating vcpkg.json..."

$jsonContent = Get-Content $vcpkgJson -Raw
$jsonContent = $jsonContent -replace '"version-date":\s*"[^"]*"', "`"version-date`": `"$versionDate`""
$jsonContent = $jsonContent -replace ',\s*"port-version":\s*\d+', ''
$jsonContent = $jsonContent -replace '"port-version":\s*\d+,?\s*', ''
Set-Content -Path $vcpkgJson -Value $jsonContent -NoNewline

Write-Host "  version-date set to $versionDate"

# --- Update portfile.cmake tag ---
Write-Host "`nUpdating portfile.cmake tag..."

$portContent = Get-Content $portfile -Raw
$portContent = $portContent -replace '(\n\s+REF\s+)\S+', "`${1}$Tag"
Set-Content -Path $portfile -Value $portContent -NoNewline

Write-Host "  Tag set to $Tag"

# --- Download and hash source archive ---
$ProgressPreference = 'SilentlyContinue'
$tempDir = Join-Path $Env:Temp $(New-Guid)
New-Item -Type Directory -Path $tempDir | Out-Null

$sourceUrl = "https://github.com/Microsoft/DirectXMath/archive/refs/tags/$Tag.tar.gz"
$sourcePath = Join-Path $tempDir "$Tag.tar.gz"

Write-Host "`nDownloading source archive from $sourceUrl..."
try {
    Invoke-WebRequest -Uri $sourceUrl -OutFile $sourcePath -ErrorAction Stop
}
catch {
    Write-Error "ERROR: Failed to download source archive!" -ErrorAction Stop
}

$sourceHash = (Get-FileHash -Path $sourcePath -Algorithm SHA512).Hash.ToLower()
Write-Host "  Source SHA512: $sourceHash"

# Replace SHA512 in vcpkg_from_github block
$portContent = Get-Content $portfile -Raw
$portContent = $portContent -replace '(vcpkg_from_github\s*\([^)]*SHA512\s+)[0-9a-fA-F]+', "`${1}$sourceHash"
Set-Content -Path $portfile -Value $portContent -NoNewline

# --- Cleanup ---
Remove-Item -Recurse -Force $tempDir

Write-Host "`nvcpkg port updated successfully!"
Write-Host "`nUpdated files:"
Write-Host "  $vcpkgJson"
Write-Host "  $portfile"

$portContent = Get-Content $portfile -Raw
if ($portContent -match '\bPATCHES\b') {
    Write-Warning "This port includes patches. Review them to either remove or update."
}
