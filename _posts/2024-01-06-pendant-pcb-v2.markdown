---
title:  "Pendant PCB v2"
date:   2024-01-06 0:00:00 +0200
categories: cnc
tags: cnc pendant
---

TLDR; Bricked a 2nd Pico for this project.

Designed a simple single layer PCB for the display and the Pico.

![PCB](/assets/img/pendant_pcb.png)

Milled and soldered. Everything worked correctly. For 5 mins.
After that, the USB device kept getting disconnected and reconnected.

After a lot of experimentation of what might went wrong, it seems that
the Pico was bricked when desoldering the pin headers. To rule out any
possible bugs in my code I loaded the USB device_hid_composite example
from the SDK. The device behaves the same.

Ordered another Pico as well as a [Waveshare RP2040 Tiny](https://www.waveshare.com/rp2040-tiny.htm) board.
The RP2040 Tiny might be a better choice for a v3 of the PCB as it allows for arbitrary placement
of the USB port on the pendant's enclosure.
