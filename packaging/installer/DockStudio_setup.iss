; Inno Setup script for DockStudio (Eric Studio)
; Build:   build_installer_windows.bat
; Output:  dist\installer\DockStudio_Setup_1.1.0.exe
;
; Requires:
;   - Inno Setup 6 (https://jrsoftware.org/isinfo.php)  -> ISCC.exe on PATH
;   - conda-pack portable folder at {#PORT} (created by build_win_condapack.bat)

#ifndef PORT
  #define PORT "..\..\dist\DockStudioPortable"
#endif
#ifndef APP_VERSION
  #define APP_VERSION "1.1.0"
#endif

#define MyAppName "DockStudio"
#define MyAppPublisher "Eric Studio"
#define MyAppURL "https://github.com/EricStudio/DockStudio"
#define MyAppExeName "DockStudio.bat"
#define MyAppCopyright "Copyright (c) 2026 Eric Studio. All rights reserved."

[Setup]
AppId={{8E3F6C8A-2C74-4E4E-B1E3-8E9D0A2D6F42}}
AppName={#MyAppName}
AppVersion={#APP_VERSION}
AppVerName={#MyAppName} {#APP_VERSION}
AppPublisher={#MyAppPublisher}
AppCopyright={#MyAppCopyright}
VersionInfoCompany=Eric Studio
VersionInfoDescription=DockStudio - one-stop batch molecular docking platform
VersionInfoTextVersion={#APP_VERSION}
VersionInfoCopyright={#MyAppCopyright}
DefaultDirName={autopf}\DockStudio
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
LicenseFile=license.txt
OutputDir=..\..\dist\installer
OutputBaseFilename=DockStudio_Setup_{#APP_VERSION}
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
SetupIconFile=..\..\assets\icon.ico
UninstallDisplayIcon={app}\dockstudio\resources\icon.ico
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64

[Languages]
Name: "chinesesimplified"; MessagesFile: "compiler:Languages\ChineseSimplified.isl"
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "{#PORT}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Comment: "DockStudio by Eric Studio"; IconFilename: "{app}\dockstudio\resources\icon.ico"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon; IconFilename: "{app}\dockstudio\resources\icon.ico"

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(MyAppName, '&', '&&')}}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
Type: filesandordirs; Name: "{app}\__pycache__"
Type: filesandordirs; Name: "{app}\dockstudio\__pycache__"
