# GAGfinder

Welcome to GAGfinder! This software was developed in Joseph Zaia's lab in the Center for Biomedical Mass Spectrometry at the Boston University School of Medicine. GAGfinder analyzes glycosaminoglycan (GAG) tandem mass spectra by searching the spectra for fragments of the determined precursor composition. Currently, GAGfinder is only available in command line form, but a web application is under development. Continue reading this file for instructions on using GAGfinder.

# Motivation

Current methods for peak finding in GAG tandem mass spectrometry rely on the use of averagine, or the average building block of the structure. Due to the variable sulfation patterns along GAG chains and the relatively strong sulfur-34 peak, these methods have a high rate of false positives and false negatives when employed on GAG tandem mass spectra. We decided to invert the steps in traditional peak finding methods and developed a targeted peak finding algorithm. GAGfinder returns GAG-related fragments at various charge states according to their determined precursor composition.

# Installation

If you've made it this far, then you have already installed GAGfinder. Congratulations!

# Running GAGfinder command line interface

In order to run GAGfinder, you need to open a Command Prompt and make your way to this current directory (gagfinder-cli_windows). If you are new to the Command Prompt, you can read a tutorial about different commands at (http://www.cs.princeton.edu/courses/archive/spr05/cos126/cmd-prompt.html). For instance, if you left the gagfinder-cli_windows folder in your Downloads folder, you could execute the following commands to get into the proper directory and list the files and folders in this folder (everything after the '>' symbol):

[prompt]>C:<br />
[prompt]>cd Users\[your username]\Downloads\gagfinder-cli_windows<br />
[prompt]>dir

Now you are in this directory. However, the GAGfinder executable is in a sub-directory of this directory, called gagfinder-cli. You will need to execute the following command to get to that directory and list the files and folders in it (everything after the '>' symbol):

[prompt]>cd gagfinder-cli<br />
[prompt]>dir

Amongst all of these files and directories should be an executable file called 'gagfinder-cli.exe'. This is your executable file to run GAGfinder. Execute the following command to ensure that the installation was successful.

[prompt]>gagfinder-cli.exe --help

You should now see in your Command Prompt a list of arguments. These are important inputs from you that tell GAGfinder exactly how to run.

# GAGfinder arguments

Below are descriptions of all of the arguments for GAGfinder:

-h, --help    This argument brings up the help screen describing all of the arguments for GAGfinder<br />
-c            This required argument is for the GAG class of your spectrum. Either HS, CS, or KS are acceptable answers. If your GAG
              class is dermatan sulfate, use CS; if your GAG class is heparin or hyaluronan, use HS.<br />
-i            This required argument is for the input mzML file. Make sure you have converted your datafile into mzML format before
              input.<br />
-r            This optional argument is for reducing end derivatization. If you have added a reducing end tag to your GAG, input its
              formula here. PLEASE NOTE: The proper way to encode your reducing end tag is to enter the chemical formula that is ADDED
              ON to the GAG. For instance, in Arixtra, the reducing end tag is listed as OMe/OCH3. However, the difference between that
              structure non-tagged and tagged is only CH2; therefore, CH2 is the proper reducing end tag to input to GAGfinder.<br />
-n            This optional argument is for the number of top results to return. Either this or the next argument (-p for top
              percentile) is required. If you are interested in receiving the top 20 results, use -n 20 as your argument.<br />
-p            This optional argument is for the percent of top results to return. Either this or the previous argument (-n for top
              number) is required. If you are interested in receiving the top 10% of results, use -p 10 as your argument.<br />
-a            This optionalargument is for the metal that is adducted. If you have adducted metal cations to your structure, you will
              need to use this argument. Only Na, K, Li, Ca, and Mg are accepted as metal adducts.<br />
-t            This optional argument is for the number of adducts added. If you input an argument for metal adduct, you must input the
              number.<br />
-g            This optional argument is for the cation reagent used in negative electron transfer dissociation. If you are interested
              in the peaks formed by the precursor binding to this reagent, enter its chemical formula.<br />
-m            This required argument is for the precursor m/z value.<br />
-z            This required argument is for the precursor charge. If the precursor charge is -3, use -z -3 as your argument.<br />
-s            This optional argument is for the number of sulfate losses to consider.<br />
-x            This required argument is for whether the noise has been removed prior to GAGfinder analysis or not. Enter y for yes or
              n for no.

# Example command

Now that you know where the executable file is and know about all of the arguments, let's go through an example analysis. Let's say you are analyzing an HS spectrum with precursor m/z of 362.4143 and precursor charge of -5. You ran the mass spectrometer in NETD mode, using the cation reagent fluoranthene (C16H10). The mzML file is located at C:\Example\Not\Real\test.mzML. The reducing end tag adds an ethyl group (-CH2CH3). You are interested in the top 20% of results, and believe that there may be up to two sulfate losses. The noise has not been removed. You would use the following command (everything after the '>' symbol):

[prompt]>gagfinder-cli.exe -c HS -i C:\Example\Not\Real\test.mzML -p 20 -g C16H10 -m 362.4143 -z -5 -s 2 -x n

The output would then be found at C:\Example\Not\Real\test.txt.

# Output format

The output includes columns for m/z, Intensity, Charge, Fragments, G-score, and Error (ppm). The G test is how GAGfinder tests each theoretical fragment's fit with what exists in the spectrum, and a lower G-score represents a better fit. The output file is sorted by ascending G-score. The fragments column shows the fragments associated with the given m/z and charge, and they have a unique string identifier. This string identifier has the following representations:

D: delta-4,5-unsaturated hexuronic acid<br />
U: Uronic acid (GlcA or IdoA)<br />
X: Hexose (Gal)<br />
N: Hexosamine (GalN or GlcN)<br />
A: Acetyl group<br />
S: Sulfate

For example, if your fragment says 'UNS3', it is a uronic acid attached to a hexosamine with three sulfate groups attached. Water loss and hydrogen loss are considered as well, so if you see a fragment that says 'UNS-H2O-H', it is a uronic acid attached to a hexosamine with one sulfate group attached, missing a water and a hydrogen. Derivative fragments of the precursor are listed as 'M', so if your fragment says 'M-SO3', it is the precursor minus one sulfate. Cross-ring cleavages are added to the string with this form:'+[monosaccharide][NR/RE][cleavage]'. For example, if your fragment says NS+UNR3,5, then it is a non-reducing end fragment containing a full hexosamine that cuts into a uronic acid on its non-reducing end side at a 3,5 cleavage. This may correspond to a 3,5A2 ion, depending on sulfate loss. If there is a reducing end tag and the fragment contains it, you will see '+RE' in the string. For example, if your fragment says NS2+URE1,5+RE, then it is a reducing-end fragment containing a full hexosamine that cuts into a uronic acid on its reducing end side at a 1,5 cleavage. This may correspond to a 1,5X1 ion, depending on sulfate loss.

# Conclusion

Congratulations! You are now ready to analyze all of those GAG spectra you have been avoiding manually interpretating. Any questions can be e-mailed to jzaia@bu.edu. Happy GAG finding!
