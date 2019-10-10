#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 16:28:15 2019

@author: jorge
"""

from aip_detector import aip_detector

file_dir = '/home/jorge/Escritorio/mitdb/100'
channel_num = 0
pattern_width = 40
pattern_time_width = 120e-3
normalize = True
initial_threshold = 30
graphics = True
create_annotation = True

aip_detector(file_dir, channel_num, pattern_width, pattern_time_width, normalize, initial_threshold, graphics, create_annotation)