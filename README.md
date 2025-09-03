# LearnMathsToday Manimations

This repo contains the Manim code for the LMT Youtube channel: 

Projects that can be found inside:

basics: beginner animations to play around with -- recurring ideas in mathematics education

arc_lengths: folder for the chapter on circles, arcs, sectors, segments and their areas, as well as the use of radians.

-------------------------------------------------
## Installation

It is advised to install Manim into a virtual environment such as a conda environment. Follow the instructions from https://www.manim.community (this is the Community version, as Grant's original repo is not designed for public use). Activate the environment and then clone this repo into a location of your choice. 


### To run the animations:

Navigate to the proper folder name, eg
`
cd LearnMathsToday/arc_lengths
`

and then

`
manim -p -qh file_name.py SceneName
`


For each folder, there are multiple Python files, and each Python files have multiple classes or scenes. For instance, the scene called "SegmentChordThenArc" can be run with the command:

```
cd arc_lengths

manim -p -qh segments.py SegmentChordThenArc
```

or alternatively you can just run 

`manim -p -qh segments.py` 

and a list of all indexed scenes will pop up which have been detected in the file segments.py


The output videos can then be found inside the workspace under the folder "media".
