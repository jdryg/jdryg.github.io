---
title:  "Pendant PCB v3"
date:   2024-01-12 0:00:00 +0200
categories: cnc
tags: cnc pendant
---

Waveshare RP2040 Tiny board arrived and we made a new PCB for it. 

![Pendant PCB with Waveshare RP2040 Tiny](/assets/img/pendant_pcb_waveshare.png)

![3D view of the PCB](/assets/img/pendant_pcb_waveshare_3d.png)

I used a PIO state machine for the display's SPI (converts 565 to 888 during the DMA transfer), so it was easier to route all the pins around the board.

Also added a debug UART which was not present in the previous version. It might be handy if something
goes wrong before the CDC is initialized.

I should now focus on the software part of the project. Most of the stuff are ready. Just a couple of hal pins
to expose and finish the UI.
