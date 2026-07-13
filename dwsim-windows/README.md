[![GitHub issues](https://img.shields.io/github/issues/DanWBR/dwsim6.svg)](https://github.com/DanWBR/dwsim6/issues)
[![tickets](https://img.shields.io/badge/view-tickets-blackgray.svg)](https://sourceforge.net/p/dwsim/tickets/)
[![forums](https://img.shields.io/badge/join-the%20forums-yellowgreen.svg)](https://sourceforge.net/p/dwsim/discussion/?source=navbar)
[![wiki](https://img.shields.io/badge/visit-website-blackblue.svg)](http://dwsim.inforside.com.br)
[![donate](https://img.shields.io/badge/make%20a-donation-greenblue.svg)](https://sourceforge.net/p/dwsim/donate/)

## DWSIM - Open Source Process Simulator
Copyright 2008-2025 Daniel Wagner and contributors

DWSIM is a software for modeling, simulating, and optimizing steady-state and dynamic chemical processes.

### License

DWSIM is licensed under the GNU General Public License (GPL) Version 3.

See COPYING for more information.

### Supported Operating Systems

- Windows (64-bit x86) with .NET Framework 4.6.2 or newer
- Linux (64-bit x86) with .NET 8 Runtime or newer
- macOS 10.7 or newer

### Donations

- Patreon: https://patreon.com/dwsim
- GitHub Sponsors: https://github.com/sponsors/DanWBR
- Buy-me-a-coffee: https://www.buymeacoffee.com/dwsim
- Bitcoin tips are welcome at bc1qf37y47vfk5wzxqpyh39y7th32x6lja0h0gc383

### Compiling

- DWSIM can be compiled using Visual Studio 2019 or newer on Windows.
- To compile everything and run:
	- Open Visual Studio 2019 or 2022 and clone this repository directly from GitHub
	- Change the Build target to 'Debug/x64', 'ReleaseLinux/x64', 'ReleaseWinMac/x64' or 'ReleaseWinMac/x86'
	- Click on the Solution object and restore NuGet packages
	- Build the solution
	- Select 'DWSIM' or 'DWSIM.UI.Desktop' as the startup project
	- Run
