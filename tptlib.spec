%define _prefix /gem_base/epics/support
%define name tptlib
%define repository gemdev
%define debug_package %{nil}
%define arch %(uname -m)
%define checkout %(git log --pretty=format:'%h' -n 1) 

#These global defines are added to prevent stripping
# symbols on vxWorks cross-compiled code
# Getting 'strip' to work is probably only needed for
# building a related debug sub-package
#
# But this prevents all the strip warnings
# mrippa 20120202
%global _enable_debug_package 0
%global debug_package %{nil}
%global __os_install_post /usr/lib/rpm/brp-compress %{nil}

Summary: %{name} Package, a module for EPICS base
Name: %{name}
Version: 1.1.1
Release: 7%{?dist}
License: EPICS Open License
Group: Applications/Engineering
Source0: %{name}-%{version}.tar.gz
ExclusiveArch: %{arch}
Prefix: %{_prefix}
## You may specify dependencies here
BuildRequires: epics-base-devel = 7.0.7-0.git.f9e3717%{?dist}
BuildRequires: re2c
BuildRequires: gemini-ade
BuildRequires: slalib-devel = 1.9.7-6.git.54d124d%{?dist}
Requires: epics-base = 7.0.7-0.git.f9e3717%{?dist}
Requires: slalib = 1.9.7-6.git.54d124d%{?dist}
## Switch dependency checking off
## AutoReqProv: no

%description
This is the module %{name}.

## If you want to have a devel-package to be generated uncomment the following:
%package devel
Summary: %{name}-devel Package
Group: Development/Gemini
Requires: %{name}
%description devel
This is the module %{name}.

%prep
%setup -q 

%build
make distclean uninstall
make

%install
export DONT_STRIP=1
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT/%{_prefix}/%{name}
#cp -r dbd $RPM_BUILD_ROOT/%{_prefix}/%{name}
#cp -r bin $RPM_BUILD_ROOT/%{_prefix}/%{name}
cp -r lib $RPM_BUILD_ROOT/%{_prefix}/%{name}
cp -r include $RPM_BUILD_ROOT/%{_prefix}/%{name}
cp -r configure $RPM_BUILD_ROOT/%{_prefix}/%{name}
# find $RPM_BUILD_ROOT/%{_prefix}/%{name}/configure -name ".git" -exec rm -rf {} \;


%postun
if [ "$1" = "0" ]; then
	rm -rf %{_prefix}/%{name}
fi


%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
#   /%{_prefix}/%{name}/bin
   /%{_prefix}/%{name}/lib

%files devel
%defattr(-,root,root)
#   /%{_prefix}/%{name}/dbd
   /%{_prefix}/%{name}/include
   /%{_prefix}/%{name}/configure

%changelog
* Mon Jul 26 2021 emmakurz <emma.kurz@noirlab.edu> 1.1.1-7
- Closes issue #1
- Squash commit collectively with Tiffany, Felix and Matt. Merge ADE 2 release
  infrastructure updates from Fast Track 3.15 work into master.

* Wed Dec 30 2020 Roberto Rojas <rrojas@gemini.edu> 1.1.1-6
- fixed a typo problem

* Wed Dec 30 2020 Roberto Rojas <rrojas@gemini.edu> 1.1.1-5
- fixed a typo problem

* Wed Dec 30 2020 Roberto Rojas <rrojas@gemini.edu> 1.1.1-4
- Added slalib to the spec file

* Wed Dec 30 2020 Roberto Rojas <rrojas@gemini.edu> 1.1.1-3
- new package built with tito

* Wed Nov 18 2020 rrojas <rrojas@gemini.edu> 1.1.1
- copy based on astlib module
