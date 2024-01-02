---
title:  "Pendant UI Mockup"
date:   2024-01-02 0:00:00 +0200
categories: cnc
tags: cnc pendant
---
Mockup of the pendant UI so far.

Power tab

![Power tab mockup](/assets/img/pendant_ui_power_tab.png)

Manual Mode tab

![Manual Mode tab mockup](/assets/img/pendant_ui_manual_mode_tab.png)

I'm still not sure if Feed and Rapid overrides should be part of the Manual Mode tab.
We've never used them at this point. Usually Manual mode is for jogging and setting the work zero offsets.
It'd probably be better if they are moved to the Program tab.

The zeroing buttons next to each axis will trigger an MDI command. Unfortunately there is no way to send custom MDI commands from a LinuxCNC component to LinuxCNC
so they must be specified in the ini file. Setting the G54 work offset is done using G10 command. I should make sure the MDI command matches the HAL pin used by
the component otherwise the pendant will trigger the wrong command.

What's missing from the Manual Mode tab is the jogging step in machine units. Since we probably don't want an arbitrary step, a slider won't be appropriate for that.
A couple of buttons to increase/decrease the step to predefined values should work.