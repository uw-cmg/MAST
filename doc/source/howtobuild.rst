How to Build Sphinx Document
****************************

Here is brief introduction for sphinx documentation.
Let me start with my awful mistake. Because of that, my first impression of sphinx is horrible. 

Even virtical space is important to compile html.

The toctree directive initially is empty, and looks like this::

    .. toctree::
        :maxdepth: 2

I added title of MAST document and add one file(document) howtobuild.rst like the following::

    Welcome to MAST's documentation!
    ================================
    
    Contents:

    .. toctree::
        :maxdepth: 2
        howtobuild

After wasting quite long time, I took a look at tutorial again. I realized that there is vertical space::

    Welcome to MAST's documentation!
    ================================

    Contents:
 
    .. toctree::
        :maxdepth: 2

        howtobuild

Then you can compile it. The compile command is this::

    ~/MAST4pymatgen/doc/make html

You may already have noticed this. Indent is also very important to get correct output. So for me it's literally annoying. You can think learning how to write sphinx document is same as learning how write code. As we did in python code, for one indent level, I used four whitespaces.

As you can see, most of Sphinx Document has ''**Show Source**'' link. In other words, all Sphinx document can be good example for you. One thing I want to tell you, in order to use all functions multiple extensions are required. Don't be frustrated even though your copied document doesn't work well.
