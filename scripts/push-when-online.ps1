# 网络恢复后自动推送 main 与全部 tag。
#
# 用途：GitHub 直连被阻断时，先提交、再让这个脚本在网络恢复后自动补推。
# 只在 `github.com:443` 可达时才尝试，避免空等超时。
#
# 用法（Windows PowerShell 5.1 与 PowerShell 7 均可）：
#   powershell -NoProfile -File scripts/push-when-online.ps1
# 可选参数：
#   -IntervalSec 45   -Attempts 40   -Remote origin
param(
    [int]$IntervalSec = 45,
    [int]$Attempts = 40,
    [string]$Remote = 'origin'
)

$ErrorActionPreference = 'Continue'
$repo = Split-Path -Parent $PSScriptRoot
Set-Location $repo

function Test-GitHub {
    Test-NetConnection -ComputerName github.com -Port 443 `
        -InformationLevel Quiet -WarningAction SilentlyContinue
}

# 需要推送的 tag：本地有、远程还没有的（动态计算，不写死版本号）
function Get-UnpushedTags {
    $local = git tag --list
    $remoteTags = (git ls-remote --tags $Remote 2>$null) -join "`n"
    $out = @()
    foreach ($t in $local) {
        if ($remoteTags -notmatch [regex]::Escape("refs/tags/$t")) { $out += $t }
    }
    return $out
}

for ($i = 1; $i -le $Attempts; $i++) {
    if (-not (Test-GitHub)) {
        Write-Output ("[{0:HH:mm:ss}] try {1}/{2}: github:443 不通，等待 {3}s" -f `
            (Get-Date), $i, $Attempts, $IntervalSec)
        Start-Sleep -Seconds $IntervalSec
        continue
    }

    Write-Output ("[{0:HH:mm:ss}] try {1}: 端口可达，开始 push" -f (Get-Date), $i)
    $out = git push $Remote main 2>&1 | Out-String
    if ($out -match 'main -> main|Everything up-to-date') {
        Write-Output '  main 推送成功'
        $tags = Get-UnpushedTags
        if ($tags.Count -gt 0) {
            Write-Output ("  推送 tag: {0}" -f ($tags -join ', '))
            git push $Remote @tags 2>&1 | Out-String | Write-Output
        } else {
            Write-Output '  没有待推送的 tag'
        }
        Write-Output 'PUSH_OK'
        exit 0
    }

    Write-Output ('  main 失败: ' + (($out -split "`n" |
        Select-String 'fatal|error' | Select-Object -First 1) -join ''))
    Start-Sleep -Seconds $IntervalSec
}

Write-Output 'PUSH_FAILED: 重试次数用尽，网络仍未恢复'
exit 1
