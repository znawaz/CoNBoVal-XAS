<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html><head>
  
  <meta content="text/html; charset=ISO-8859-1" http-equiv="content-type">
  <title>CoNBoVal-XAS ReadMe</title>

  
</head><body>
<h1 style="font-family: Arial;">CoNBoVal-XAS</h1>

<span style="font-family: Arial;">CoNBoVal_XAS</span> : <span style="font-family: Arial;">Coordination Number Bond Valence for X-Ray
Absorption Spectroscopy<br>
<br>
CoNBoVal_XAS is a program which takes atom positions in three
dimensions and computes the following<br>
</span>
<ul>

  <li><span style="font-family: Arial;">Nearest neighborhood distance
of the neighboring atoms for a given atom and with in a specified
distance</span>. <span style="font-family: Arial;">This output becomes
the input for the FEFF software.</span></li>
  <li><span style="font-family: Arial;">Coordination Number when
provided with coordination atom</span></li>
  <li><span style="font-family: Arial;">Bond Valence for all the
coordination atoms.</span></li>
</ul>

<span style="font-family: Arial;">The 3D atom positions can be the
output of a DFT or molecular dynamics simulation. After computing the
nearest neighborhood distance, it becomes the input for the FEFF
software.<br>
<br>
The computation of Coordination number and Bond Valence of the
coordinating atom can be used as validation of DFT or molecular
dynamics simulation.<br>
<br>
</span>
<h2><span style="font-family: Arial;">How to Build the CoNBoVal-XAS</span></h2>

<span style="font-family: Arial;">To build the CoNBoVal-XAS, you need
to download ANN (Approximate Nearest Neighbors) library from <a href="http://www.cs.umd.edu/%7Emount/ANN">http://www.cs.umd.edu/~mount/ANN</a>.</span>&nbsp;
<span style="font-family: Arial;">Dowload the latest zip file in a
directory. Then you need to set ANN_HOME environment variable. In
linux, it can be set as follows<br>
<br>
<span style="font-family: Courier New;">export ANN_HOME=[path to ann
library]</span></span><code style="border: medium none ; margin: 0px; padding: 0px; font-family: Courier New; font-size: 12px; background-color: transparent; white-space: pre; display: inline; line-height: inherit;"></code><span style="font-family: Arial;"><span style="font-family: Courier New;">/ann_1.1.2</span><br style="font-family: Courier New;">
<br>
like on my system it is<br>
<span style="font-family: Courier New;">export
ANN_HOME=/home/zubair/software/ann_1.1.2</span><br>
<br>
If you want to keep it permanent, copy this line in .bashrc.<br>
<br>
To build, type in the root directory of CoNBoVal-XAS the following<br>
<span style="font-family: Courier New;">make</span><br>
<br>
</span>
<h2><span style="font-family: Arial;">How to use</span></h2>

<span style="font-family: Arial;">The program is invoked as follows<br>
<br><span style="font-family: Courier New;">
CoNBoVal-XAS [-i inputfile] [-d prefix] [-b bondValence]</span><br>
<br>
where:<br><span style="font-family: Courier New;">
inputfile&nbsp;</span>&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; name
of the input file containing position of atoms (default = "MD-data")<br><span style="font-family: Courier New;">
prefix&nbsp;&nbsp;</span>&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;
&nbsp;&nbsp; &nbsp; prefix name of the output files that contain
distances (default = "feff")<br><span style="font-family: Courier New;">
bondValence</span>&nbsp;&nbsp; suffix name of the output file that contain
bond valence (default = "_bondVal.txt")<br>
&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; <br>
Results are stored in prefix and bondValence files.<br>
<br><span style="font-weight: bold;">
Example:&nbsp;</span> &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;
&nbsp;&nbsp; <br>
To run this demo use:<br>
<span style="font-family: Courier New;">CoNBoVal-XAS -i MD-data -d feff -b bond-valence.txt</span><br>
or<br><span style="font-family: Courier New;">
CoNBoVal-XAS<br>
<br>
</span></span>
<h2><span style="font-family: Arial;"><span style="font-family: Courier New;"><span style="font-family: Arial;">Author</span></span></span></h2>
<span style="font-family: Arial;"><span style="font-family: Courier New;"><span style="font-family: Arial;">The main author of this program is Zubair Nawaz (<a href="mailto:zubair.nawaz@gmail.com">zubair.nawaz@gmail.com</a>). It is development with the help of Beamline Scientist Messaoud Harfouche</span>.<br>
<br>
</span></span>
<h2><span style="font-family: Arial;"><span style="font-family: Courier New;"><span style="font-family: Arial;">Acknowledgement</span></span></span></h2>
<h2><span style="font-family: Arial;"><span style="font-family: Courier New;"></span></span></h2>
<span style="font-family: Arial;"><span style="font-family: Courier New;"><span style="font-family: Arial;">&nbsp;</span><span style="font-family: Arial;">This
work was developed in SESAME and funded under the LinkSCEEM-2 project
under Grant Agreement 261600 of the 7th framework programme for
research (FP7)</span>.<br>
</span>
</span>
<h2><span style="font-family: Arial;"></span></h2>

<span style="font-family: Arial;"><br>
</span><br>

<h2><span style="font-family: Arial;"></span></h2>

<br>

<br>
</body></html>