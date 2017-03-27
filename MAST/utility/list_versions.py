#!/usr/bin/env python

#http://stackoverflow.com/questions/12939975/how-to-list-all-installed-packages-and-their-versions-in-python  frosty

import pip #needed to use the pip functions
for i in pip.get_installed_distributions(local_only=True):
    print(i)
