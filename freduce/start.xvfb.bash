#!/bin/bash

# run Xvfb to listen on ":1.0"
# to write graphics to this Xvfb instance set environment "DISPLAY=:1.0"
sudo Xvfb :1 -screen 0 1280x1024x24 &
