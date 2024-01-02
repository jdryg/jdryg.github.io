---
layout: post
title:  "Probe Basic initial tests"
date:   2024-01-01 10:00:00 +0200
categories: cnc
---
Installed Linux CNC 2.9.2 yesterday.

Rebuilt the machine in gmoccapy to make sure everything continue to work correctly. We had some issues with the Mesa card last time we tested 2.9.1.

Installed [Probe Basic][probe-basic] and followed [the guide][probe-basic-guide] on the LinuxCNC forum to configure it.
The machine seems to function correctly. 

Probing works great. That was the main reason for the transition.
The calibration tab calculated a 0.06mm deviation for our custom probe.
After calibration both axes were within 1 step (< 3.125um) away from the true value (measured the inner diameter of a ball bearing which we have no way to know the real value but we assume it's the nominal).

[probe-basic]: https://kcjengr.github.io/probe_basic/
[probe-basic-guide]: https://forum.linuxcnc.org/qtpyvcp/48401-py3-probe-basic-config-conversion-doc-lcnc-2-9
