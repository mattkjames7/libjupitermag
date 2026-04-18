param(
    [string]$BuildDir = "build-cmake"
)

$ErrorActionPreference = 'Stop'

$lmic = Join-Path $BuildDir '_deps/libcon2020-src/src/lmic.cc'
if (-not (Test-Path $lmic)) {
    throw "Missing expected file: $lmic"
}
$lmicContent = Get-Content $lmic -Raw
$lmicPatched = $lmicContent -replace 'double\s+B\s*=\s*abs\(', 'double B = fabs('
if ($lmicPatched -eq $lmicContent) {
    throw "Did not find expected abs(...) expression in: $lmic"
}
Set-Content -Path $lmic -Value $lmicPatched -NoNewline

$internal = Join-Path $BuildDir '_deps/libinternalfield-src/src/internal.cc'
if (-not (Test-Path $internal)) {
    throw "Missing expected file: $internal"
}
$internalContent = Get-Content $internal -Raw

$oldFacts = 'double facts[nfact];'
$newFacts = 'double *facts = new double[nfact];'
if (-not $internalContent.Contains($oldFacts)) {
    throw "Did not find expected VLA declaration in: $internal"
}
$internalContent = $internalContent.Replace($oldFacts, $newFacts)

$oldTail = @"
            Snm_[n][m] = sqrt(delta*((facts[n-m]/facts[n+m])));
        }
    }
}
"@
$newTail = @"
            Snm_[n][m] = sqrt(delta*((facts[n-m]/facts[n+m])));
        }
    }

    delete[] facts;
}
"@
if (-not $internalContent.Contains($oldTail)) {
    throw "Did not find expected Schmidt tail in: $internal"
}
$internalContent = $internalContent.Replace($oldTail, $newTail)

Set-Content -Path $internal -Value $internalContent -NoNewline

Write-Host "Applied MSVC compatibility patches in fetched dependencies."
