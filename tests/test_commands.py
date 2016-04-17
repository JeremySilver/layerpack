"""Tests for layerpack  CLI module."""

## forked from: https://github.com/rdegges/skele-cli

from subprocess import PIPE, Popen as popen
from unittest import TestCase

from layerpack import __version__ as VERSION

class TestHelp(TestCase):
    def test_returns_usage_information(self):
        commands = ['ncpack', 'ncunpack', 'nccheckdiff']
        helpflags = ['-h', '--help']
        ##
        for cmd in commands:
            for flg in helpflags:
                output = popen([cmd, flg], stdout=PIPE).communicate()[0]
                self.assertTrue('Usage:' in output)


class TestVersion(TestCase):
    def test_returns_version_information(self):
        commands = ['ncpack', 'ncunpack', 'nccheckdiff']
        ##
        for cmd in commands:
            output = popen([cmd, '--version'], stdout=PIPE).communicate()[0]
            self.assertEqual(output.strip(), VERSION)
                
