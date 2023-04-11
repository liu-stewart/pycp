"""This module used to load the settings."""
import json
import pathlib

with open(pathlib.Path.home()/".pycp/setting.json", "r") as file:
    setting = json.load(file)
