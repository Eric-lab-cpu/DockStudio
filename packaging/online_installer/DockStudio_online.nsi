; ============================================================
; DockStudio (Eric Studio) - Online Installer
; Build: makensis DockStudio_online.nsi   (runs on Linux or Windows)
; Output: <repo>/dist/DockStudio_Setup_1.1.0.exe
;
; The payload (dist/payload) is prepared by prepare_payload.py
; ============================================================

!include "MUI2.nsh"
!include "FileFunc.nsh"

!define APP_NAME "DockStudio"
!define APP_VERSION "1.1.0"
!define APP_PUBLISHER "Eric Studio"
!define APP_COPYRIGHT "Copyright (c) 2026 Eric Studio. All rights reserved."
!define PORT "..\..\dist\payload"
!define ICON "..\..\assets\icon.ico"
!define UNINST_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\DockStudio"

Name "${APP_NAME}"
OutFile "..\..\dist\DockStudio_Setup_${APP_VERSION}.exe"
InstallDir "$LOCALAPPDATA\EricStudio\DockStudio"
InstallDirRegKey HKCU "Software\EricStudio\DockStudio" ""
RequestExecutionLevel user
Unicode true
SetCompressor /SOLID lzma
CRCCheck on

VIProductVersion "${APP_VERSION}.0.0"
VIAddVersionKey "ProductName" "${APP_NAME}"
VIAddVersionKey "ProductVersion" "${APP_VERSION}"
VIAddVersionKey "FileVersion" "${APP_VERSION}"
VIAddVersionKey "CompanyName" "${APP_PUBLISHER}"
VIAddVersionKey "LegalCopyright" "${APP_COPYRIGHT}"
VIAddVersionKey "FileDescription" "${APP_NAME} - one-stop batch molecular docking platform"

!define MUI_ICON "${ICON}"
!define MUI_UNICON "${ICON}"
!define MUI_ABORTWARNING
!define MUI_FINISHPAGE_RUN "$INSTDIR\DockStudio.bat"
!define MUI_FINISHPAGE_RUN_TEXT "立即启动 ${APP_NAME} (首次启动会自动配置环境)"
!define MUI_FINISHPAGE_LINK "${APP_PUBLISHER}"
!define MUI_FINISHPAGE_LINK_LOCATION "https://github.com/EricStudio/DockStudio"

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "..\installer\license.txt"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES

!insertmacro MUI_LANGUAGE "SimpChinese"
!insertmacro MUI_LANGUAGE "English"
!insertmacro MUI_RESERVEFILE_LANGDLL

Section "DockStudio (required)" SecMain
    SectionIn RO
    SetOutPath "$INSTDIR"
    File /r "${PORT}\*.*"

    WriteUninstaller "$INSTDIR\Uninstall.exe"
    WriteRegStr HKCU "Software\EricStudio\DockStudio" "" "$INSTDIR"
    WriteRegStr HKCU "${UNINST_KEY}" "DisplayName" "${APP_NAME} ${APP_VERSION}"
    WriteRegStr HKCU "${UNINST_KEY}" "DisplayVersion" "${APP_VERSION}"
    WriteRegStr HKCU "${UNINST_KEY}" "Publisher" "${APP_PUBLISHER}"
    WriteRegStr HKCU "${UNINST_KEY}" "DisplayIcon" "$INSTDIR\icon.ico"
    WriteRegStr HKCU "${UNINST_KEY}" "UninstallString" '"$INSTDIR\Uninstall.exe"'
    WriteRegStr HKCU "${UNINST_KEY}" "InstallLocation" "$INSTDIR"
    WriteRegStr HKCU "${UNINST_KEY}" "NoModify" "1"
    WriteRegStr HKCU "${UNINST_KEY}" "NoRepair" "1"
    WriteRegDWORD HKCU "${UNINST_KEY}" "EstimatedSize" 65536

    CreateDirectory "$SMPROGRAMS\DockStudio"
    CreateShortCut "$SMPROGRAMS\DockStudio\DockStudio.lnk" "$INSTDIR\DockStudio.bat" "" "$INSTDIR\icon.ico"
    CreateShortCut "$SMPROGRAMS\DockStudio\卸载 DockStudio.lnk" "$INSTDIR\Uninstall.exe" "" "$INSTDIR\icon.ico"
    CreateShortCut "$DESKTOP\DockStudio.lnk" "$INSTDIR\DockStudio.bat" "" "$INSTDIR\icon.ico"
SectionEnd

Section "Uninstall"
    Delete "$INSTDIR\Uninstall.exe"
    RMDir /r "$INSTDIR\env"
    RMDir /r "$INSTDIR\.pkgs"
    RMDir /r "$INSTDIR\.mamba"
    RMDir /r "$INSTDIR"
    Delete "$SMPROGRAMS\DockStudio\DockStudio.lnk"
    Delete "$SMPROGRAMS\DockStudio\卸载 DockStudio.lnk"
    RMDir "$SMPROGRAMS\DockStudio"
    Delete "$DESKTOP\DockStudio.lnk"
    DeleteRegKey HKCU "${UNINST_KEY}"
    DeleteRegKey HKCU "Software\EricStudio\DockStudio"
SectionEnd
