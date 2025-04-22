# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2022 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import sys      # for argv[] if called as separate process / program
import os       # for os.path
from tkinter import Toplevel, BooleanVar, W, Checkbutton, Label, LabelFrame, Frame
try:
    from tkmacosx import Button
except ImportError:
    from tkinter import Button

from utilities import nativeOS, DEBUG, unquote, wgse_message
from commandprocessor import run_bash_script
# from mainwindow import mainwindow_resume      # Localized inside cancel_autosomal_formats_window() due to loop
import settings as wgse


# Todo Microarray Format Window update -- as tab to main window? New standalone program and class?

"""###################################################################################################################
  Microarray File Export subsystem / window (aka Autosomal File Formats)
    This microarray files generator developed with WGSE Extract by City Farmer.  The largest, most
    complete and unique contribution of the WGS Extract tool.  Grew out of the basic concept in Extract23 that took a
    30x WGS and generated a 23andMe v3 file.  Now based on a generated CombinedKit template file (in 37/38 reference
    format) which is then subsetted for each of 12+ formats (using a new module "aconv" developed here). Also does
    a liftover (based on wrapping pyliftover rewrite of UCSC's liftover) from Build 38 to 37.
"""
checkbuttons_results = []

# IMPORTANT: Names must match exactly as used in Button label definition ; order is not important
kits = [
    'CombinedKit', '23andMe_V3', '23andMe_V4', '23andMe_V5', '23andMe_SNPs_API', '23andMe_V35',
    'Ancestry_V1', 'Ancestry_V2', 'FTDNA_V2', 'FTDNA_V3', 'LDNA_V1', 'LDNA_V2',
    'MyHeritage_V1', 'MyHeritage_V2', "MTHFRGen", "Genera", "meuDNA", "1240K", "HOv1", "1240+HO"
]

# Any use of these individual variables is not order important. Just kept the same order for consistency
combinedButton      = None   ;   combinedButtonResult      = None
g23andmev3Button    = None   ;   g23andmev3ButtonResult    = None
g23andmev4Button    = None   ;   g23andmev4ButtonResult    = None
g23andmev5Button    = None   ;   g23andmev5ButtonResult    = None
g23andmeAPIButton   = None   ;   g23andmeAPIButtonResult   = None
g23andmev35Button   = None   ;   g23andmev35ButtonResult   = None
ancestryv1Button    = None   ;   ancestryv1ButtonResult    = None
ancestryv2Button    = None   ;   ancestryv2ButtonResult    = None
ftdnav2Button       = None   ;   ftdnav2ButtonResult       = None
ftdnav3Button       = None   ;   ftdnav3ButtonResult       = None
ldnav1Button        = None   ;   ldnav1ButtonResult        = None
ldnav2Button        = None   ;   ldnav2ButtonResult        = None
myheritagev1Button  = None   ;   myheritagev1ButtonResult  = None
myheritagev2Button  = None   ;   myheritagev2ButtonResult  = None
mthfrGenButton      = None   ;   mthfrGenButtonResult      = None
generaButton        = None   ;   generaButtonResult        = None
meuDNAButton        = None   ;   meuDNAButtonResult        = None
g1240KButton        = None   ;   g1240KButtonResult        = None
hov1Button          = None   ;   hov1ButtonResult          = None
g1240hoButton       = None   ;   g1240hoButtonResult       = None

selectAutosomalFormatsWindow = None
generateSelectedFilesButton = None

minMicZipsize = 2500000         # Minimum generated Microarray ZIP size; pretty stable at size of template
minCbnKitSize = 500000          # Parameterized minimum Combined Kit size (zip'ped) that is considered valid (Bytes)
# Most often is 17.5 MB. Had it set at 5 MB but Teemu wants it smaller for aDNA samples he runs through to get some
#  ad-mixture analysis from. Cannot recall why was set at 5 MB; maybe due to some common error we were encountering.


def button_select_autosomal_formats():
    """
        Main entry into the microarray file generator from WGSE (or as a standalone with locally recreated settings).
        Generates main selection window and then waits on User direction. Utilizes selected (in settings) BAM and
        Output directory area.  Final output are microarray files selected by user.
        Note: if finds an existing CombinedKit file, will start with that instead of recreating it
    """
    global combinedButton
    global g23andmev3Button, g23andmev4Button, g23andmev5Button, g23andmeAPIButton, g23andmev35Button
    global ancestryv1Button, ancestryv2Button, ftdnav2Button, ftdnav3Button
    global ldnav1Button, ldnav2Button, myheritagev1Button, myheritagev2Button
    global mthfrGenButton, generaButton, meuDNAButton, g1240KButton, hov1Button, g1240hoButton
    global combinedButtonResult
    global g23andmev3ButtonResult, g23andmev4ButtonResult, g23andmev5ButtonResult, g23andmeAPIButtonResult, g23andmev35ButtonResult
    global ancestryv1ButtonResult, ancestryv2ButtonResult, ftdnav2ButtonResult, ftdnav3ButtonResult
    global ldnav1ButtonResult, ldnav2ButtonResult, myheritagev1ButtonResult, myheritagev2ButtonResult
    global mthfrGenButtonResult, generaButtonResult, meuDNAButtonResult, g1240KButtonResult, hov1ButtonResult, g1240hoButtonResult
    global generateSelectedFilesButton, checkbuttons_results, selectAutosomalFormatsWindow

    font = wgse.fonts.table

    # Create main selection window for the Autosomal microarray file generator
    # wgse.window.withdraw()
    selectAutosomalFormatsWindow = Toplevel(wgse.window)
    selectAutosomalFormatsWindow.transient()
    selectAutosomalFormatsWindow.protocol("WM_DELETE_WINDOW", cancel_autosomal_formats_window)
    selectAutosomalFormatsWindow.title(wgse.lang.i18n['TitleSelectAutosomalFormats'])
    selectAutosomalFormatsWindow.geometry(wgse.microarr_winsize)
    selectAutosomalFormatsWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    selectAutosomalFormatsWindow.columnconfigure(0, weight=1)
    selectAutosomalFormatsWindow.rowconfigure(0, weight=1)

    questionAutosomalFormats = Label(selectAutosomalFormatsWindow, text=wgse.lang.i18n['QuestionWhichAutosomalFormats'],
                                     font=font['16b'], anchor="nw", justify="left")
    questionAutosomalFormats.grid(column=0, row=0, padx=5, pady=5)

    # Create Frames for each company format and row of buttons across the bottom
    combinedFrame = LabelFrame(selectAutosomalFormatsWindow, text=wgse.lang.i18n['FrameEverything'], font=font['14'])
    combinedFrame.grid(row=1, padx=5, pady=3, sticky='W')
    combinedFrame.columnconfigure(0, weight=1)
    combinedFrame.rowconfigure(1, weight=1)
    g23andmeFrame = LabelFrame(selectAutosomalFormatsWindow, text=wgse.lang.i18n['Frame23andMe'], font=font['14'])
    g23andmeFrame.grid(row=2, padx=5, pady=3, sticky='W')
    g23andmeFrame.columnconfigure(0, weight=1)
    g23andmeFrame.rowconfigure(1, weight=1)
    ancestryFrame = LabelFrame(selectAutosomalFormatsWindow, text=wgse.lang.i18n['FrameAncestry'], font=font['14'])
    ancestryFrame.grid(row=3, padx=5, pady=3, sticky='W')
    ancestryFrame.columnconfigure(0, weight=1)
    ancestryFrame.rowconfigure(1, weight=1)
    ftdnaFrame = LabelFrame(selectAutosomalFormatsWindow, text=wgse.lang.i18n['FrameFamilyTreeDNA'], font=font['14'])
    ftdnaFrame.grid(row=4, padx=5, pady=3, sticky='W')
    ftdnaFrame.columnconfigure(0, weight=1)
    ftdnaFrame.rowconfigure(1, weight=1)
    ldnaFrame = LabelFrame(selectAutosomalFormatsWindow, text=wgse.lang.i18n['FrameLivingDNA'], font=font['14'])
    ldnaFrame.grid(row=5, padx=5, pady=3, sticky='W')
    ldnaFrame.columnconfigure(0, weight=1)
    ldnaFrame.rowconfigure(1, weight=1)
    myheritageFrame = LabelFrame(selectAutosomalFormatsWindow, text=wgse.lang.i18n['FrameMyHeritage'], font=font['14'])
    myheritageFrame.grid(row=6, padx=5, pady=3, sticky='W')
    myheritageFrame.columnconfigure(0, weight=1)
    myheritageFrame.rowconfigure(1, weight=1)
    othersFrame = LabelFrame(selectAutosomalFormatsWindow, text=wgse.lang.i18n['FrameOthers'], font=font['14'])
    othersFrame.grid(row=7, padx=5, pady=3, sticky='W') if wgse.DEBUG_MODE else othersFrame.grid_remove()
    othersFrame.columnconfigure(0, weight=1)
    othersFrame.rowconfigure(1, weight=1)
    reichFrame = LabelFrame(selectAutosomalFormatsWindow, text=wgse.lang.i18n['FrameReichLab'], font=font['14'])
    reichFrame.grid(row=8, padx=5, pady=3, sticky='W') if wgse.DEBUG_MODE else othersFrame.grid_remove()
    reichFrame.columnconfigure(0, weight=1)
    reichFrame.rowconfigure(1, weight=1)
    buttonsFrame = Frame(selectAutosomalFormatsWindow)
    buttonsFrame.grid(row=9, padx=5, pady=3, sticky='W')
    buttonsFrame.columnconfigure(0, weight=1)
    buttonsFrame.rowconfigure(1, weight=1)

    # Others are MTHFR Genetics (UK), Genera (Brazil), meuDNA (Brazil)
    # Reich (Lab, Allen Ancient DNA Resource) are 1240K and Human Origins (HO)

    # Generate select box buttons for each available format
    checkbuttons_results = []

    combinedSpacer = Label(combinedFrame, text="     ", font=font['14'], anchor="w")
    combinedSpacer.grid(row=0, column=0, padx=5, pady=2)

    combinedButtonResult = BooleanVar()
    combinedButtonResult.set(False)
    checkbuttons_results.append(False)
    combinedButton = Checkbutton(combinedFrame, variable=combinedButtonResult, font=font['14'], anchor=W,
                text=wgse.lang.i18n['CombinedFileAllSNPs'],
                command=lambda: adna_checkbutton_click("CombinedKit"))
    combinedButton.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    combinedTrailingSpacer = Label(combinedFrame, text="     ", font=font['14'], anchor="w")
    combinedTrailingSpacer.grid(row=0, column=2, padx=5, pady=2)

    g23andmeSpacer = Label(g23andmeFrame, text="     ", font=font['14'], anchor="w")
    g23andmeSpacer.grid(row=0, column=0, padx=5, pady=2)

    g23andmev3ButtonResult = BooleanVar()
    g23andmev3ButtonResult.set(False)
    checkbuttons_results.append(False)
    g23andmev3Button = Checkbutton(g23andmeFrame, variable=g23andmev3ButtonResult, font=font['14'], anchor=W,
                text='v3 (11/2010 - 11/2013) (' + wgse.lang.i18n['RECOMMENDED'] + ')',
                command=lambda: adna_checkbutton_click("23andMe_V3"))
    g23andmev3Button.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    g23andmev4ButtonResult = BooleanVar()
    g23andmev4ButtonResult.set(False)
    checkbuttons_results.append(False)
    g23andmev4Button = Checkbutton(g23andmeFrame, variable=g23andmev4ButtonResult, font=font['14'], anchor=W,
                text='v4 (11/2013 - 08/2017)',
                command=lambda: adna_checkbutton_click("23andMe_V4"))
    g23andmev4Button.grid(row=0, column=2, padx=1, pady=1, sticky='W')

    g23andmev5ButtonResult = BooleanVar()
    g23andmev5ButtonResult.set(False)
    checkbuttons_results.append(False)
    g23andmev5Button = Checkbutton(g23andmeFrame, variable=g23andmev5ButtonResult, font=font['14'], anchor=W,
                text='v5 (' + wgse.lang.i18n['Since'] + ' 08/2017) (' + wgse.lang.i18n['RECOMMENDED'] + ')',
                command=lambda: adna_checkbutton_click("23andMe_V5"))
    g23andmev5Button.grid(row=1, column=1, padx=1, pady=1, sticky='W')

    g23andmeAPIButtonResult = BooleanVar()
    g23andmeAPIButtonResult.set(False)
    checkbuttons_results.append(False)
    g23andmeAPIButton = Checkbutton(g23andmeFrame, variable=g23andmeAPIButtonResult, font=font['14'], anchor=W,
                text=wgse.lang.i18n['23andMeFutureSNPs'],
                command=lambda: adna_checkbutton_click("23andMe_SNPs_API"))
    g23andmeAPIButton.grid(row=1, column=2, padx=1, pady=1, sticky='W')

    g23andmev35ButtonResult = BooleanVar()
    g23andmev35ButtonResult.set(False)
    checkbuttons_results.append(False)
    g23andmev35Button = Checkbutton(g23andmeFrame, variable=g23andmev35ButtonResult, font=font['14'], anchor=W,
                text='v3 + v5 (' + wgse.lang.i18n['RECOMMENDED'] + ')',
                command=lambda: adna_checkbutton_click("23andMe_V35"))
    g23andmev35Button.grid(row=2, column=1, padx=1, pady=1, sticky='W')

    ancestrySpacer = Label(ancestryFrame, text="     ", font=font['14'], anchor="w")
    ancestrySpacer.grid(row=0, column=0, padx=5, pady=2)

    ancestryv1ButtonResult = BooleanVar()
    ancestryv1ButtonResult.set(False)
    checkbuttons_results.append(False)
    ancestryv1Button = Checkbutton(ancestryFrame, variable=ancestryv1ButtonResult, font=font['14'], anchor=W,
                text='v1 (01/2012 - 05/2016)',
                command=lambda: adna_checkbutton_click("Ancestry_V1"))
    ancestryv1Button.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    ancestryv2ButtonResult = BooleanVar()
    ancestryv2ButtonResult.set(False)
    checkbuttons_results.append(False)
    ancestryv2Button = Checkbutton(ancestryFrame, variable=ancestryv2ButtonResult, font=font['14'], anchor=W,
                text='v2 (' + wgse.lang.i18n['Since'] + ' 05/2016)',
                command=lambda: adna_checkbutton_click("Ancestry_V2"))
    ancestryv2Button.grid(row=0, column=2, padx=1, pady=1, sticky='W')

    ftdnaSpacer = Label(ftdnaFrame, text="     ", font=font['14'], anchor="w")
    ftdnaSpacer.grid(row=0, column=0, padx=5, pady=2)

    ftdnav2ButtonResult = BooleanVar()
    ftdnav2ButtonResult.set(False)
    checkbuttons_results.append(False)
    ftdnav2Button = Checkbutton(ftdnaFrame, variable=ftdnav2ButtonResult, font=font['14'], anchor=W,
                text='v2 (02/2011 - 04/2019)',
                command=lambda: adna_checkbutton_click("FTDNA_V2"))
    ftdnav2Button.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    ftdnav3ButtonResult = BooleanVar()
    ftdnav3ButtonResult.set(False)
    checkbuttons_results.append(False)
    ftdnav3Button = Checkbutton(ftdnaFrame, variable=ftdnav3ButtonResult, font=font['14'], anchor=W,
                text='v3 (' + wgse.lang.i18n['Since'] + ' 04/2019)',
                command=lambda: adna_checkbutton_click("FTDNA_V3"))
    ftdnav3Button.grid(row=0, column=2, padx=1, pady=1, sticky='W')

    ldnaSpacer = Label(ldnaFrame, text="     ", font=font['14'], anchor="w")
    ldnaSpacer.grid(row=0, column=0, padx=5, pady=2)

    ldnav1ButtonResult = BooleanVar()
    ldnav1ButtonResult.set(False)
    checkbuttons_results.append(False)
    ldnav1Button = Checkbutton(ldnaFrame, variable=ldnav1ButtonResult, font=font['14'], anchor=W,
                text='v1 (09/2016 - 10/2018)',
                command=lambda: adna_checkbutton_click("LDNA_V1"))
    ldnav1Button.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    ldnav2ButtonResult = BooleanVar()
    ldnav2ButtonResult.set(False)
    checkbuttons_results.append(False)
    ldnav2Button = Checkbutton(ldnaFrame, variable=ldnav2ButtonResult, font=font['14'], anchor=W,
                text='v2 (' + wgse.lang.i18n['Since'] + ' 10/2018)',
                command=lambda: adna_checkbutton_click("LDNA_V2"))
    ldnav2Button.grid(row=0, column=2, padx=1, pady=1, sticky='W')

    myheritageSpacer = Label(myheritageFrame, text="     ", font=font['14'], anchor="w")
    myheritageSpacer.grid(row=0, column=0, padx=5, pady=2)

    myheritagev1ButtonResult = BooleanVar()
    myheritagev1ButtonResult.set(False)
    checkbuttons_results.append(False)
    myheritagev1Button = Checkbutton(myheritageFrame, variable=myheritagev1ButtonResult, font=font['14'], anchor=W,
                text='v1 (11/2016 - 03/2019)',
                command=lambda: adna_checkbutton_click("MyHeritage_V1"))
    myheritagev1Button.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    myheritagev2ButtonResult = BooleanVar()
    myheritagev2ButtonResult.set(False)
    checkbuttons_results.append(False)
    myheritagev2Button = Checkbutton(myheritageFrame, variable=myheritagev2ButtonResult, font=font['14'], anchor=W,
                text='v2 (' + wgse.lang.i18n['Since'] + ' 03/2019)',
                command=lambda: adna_checkbutton_click("MyHeritage_V2"))
    myheritagev2Button.grid(row=0, column=2, padx=1, pady=1, sticky='W')

    othersSpacer = Label(othersFrame, text="     ", font=font['14'], anchor="w")
    othersSpacer.grid(row=0, column=0, padx=5, pady=2)

    mthfrGenButtonResult = BooleanVar()
    mthfrGenButtonResult.set(False)
    checkbuttons_results.append(False)
    mthfrGenButton = Checkbutton(othersFrame, variable=mthfrGenButtonResult, font=font['14'], anchor=W,
                text='MTHFR Genetics (UK)',
                command=lambda: adna_checkbutton_click("MTHFRGen"))
    mthfrGenButton.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    generaButtonResult = BooleanVar()
    generaButtonResult.set(False)
    checkbuttons_results.append(False)
    generaButton = Checkbutton(othersFrame, variable=generaButtonResult, font=font['14'], anchor=W,
                text='Genera (BR)',
                command=lambda: adna_checkbutton_click("Genera"))
    generaButton.grid(row=1, column=1, padx=1, pady=1, sticky='W')

    meuDNAButtonResult = BooleanVar()
    meuDNAButtonResult.set(False)
    checkbuttons_results.append(False)
    meuDNAButton = Checkbutton(othersFrame, variable=meuDNAButtonResult, font=font['14'], anchor=W,
                text='meuDNA (BR)',
                command=lambda: adna_checkbutton_click("meuDNA"))
    meuDNAButton.grid(row=1, column=2, padx=1, pady=1, sticky='W')

    reichSpacer = Label(reichFrame, text="     ", font=font['14'], anchor="w")
    reichSpacer.grid(row=0, column=0, padx=5, pady=2)

    g1240KButtonResult = BooleanVar()
    g1240KButtonResult.set(False)
    checkbuttons_results.append(False)
    g1240KButton = Checkbutton(reichFrame, variable=g1240KButtonResult, font=font['14'], anchor=W,
                text='AADR 1240K',
                command=lambda: adna_checkbutton_click("1240K"))
    g1240KButton.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    hov1ButtonResult = BooleanVar()
    hov1ButtonResult.set(False)
    checkbuttons_results.append(False)
    hov1Button = Checkbutton(reichFrame, variable=hov1ButtonResult, font=font['14'], anchor=W,
                text='Human Origins v1',
                command=lambda: adna_checkbutton_click("HOv1"))
    hov1Button.grid(row=0, column=2, padx=1, pady=1, sticky='W')

    g1240hoButtonResult = BooleanVar()
    g1240hoButtonResult.set(False)
    checkbuttons_results.append(False)
    g1240hoButton = Checkbutton(reichFrame, variable=g1240hoButtonResult, font=font['14'], anchor=W,
                text='1240K + Human Origins v1',
                command=lambda: adna_checkbutton_click("1240+HO"))
    g1240hoButton.grid(row=1, column=1, padx=1, pady=1, sticky='W')

    # Generate action buttons across the bottom row
    selectAllFileFormatsButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonSelectEverything'],
                font=font['14'], command=button_select_every_autosomal_test)
    selectAllFileFormatsButton.grid(column=0, row=0, padx=5, pady=3)

    selectRecFileFormatsButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonSelectRecommended'],
                font=font['14'], command=button_select_recommended_autosomal_test)
    selectRecFileFormatsButton.grid(column=1, row=0, padx=5, pady=3)

    deSelectAllFileFormatsButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonDeselectEverything'],
                font=font['14'], command=button_deselect_every_autosomal_test)
    deSelectAllFileFormatsButton.grid(column=2, row=0, padx=5, pady=3)

    generateSelectedFilesButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonGenerateSelectedFiles'],
                font=font['14'], command=button_generate_selected_autosomal, state="disabled")
    generateSelectedFilesButton.grid(column=4, row=0, padx=5, pady=3)

    # Had previously assumed user knew they could click window exit in upper bar; now provide explicit exit button
    exitButton = Button(buttonsFrame, text=wgse.lang.i18n['CloseWindow'], font=font['14'],
                        command=cancel_autosomal_formats_window)
    exitButton.grid(column=5, row=0, padx=5, pady=3)

    # Select Recommended by default and wait for an action button (generate or exit)
    button_select_recommended_autosomal_test()
    selectAutosomalFormatsWindow.update()
    selectAutosomalFormatsWindow.deiconify()    # Update not showing recommended at start after grab_set() added
    selectAutosomalFormatsWindow.grab_set()
    selectAutosomalFormatsWindow.wait_window()


# noinspection PyUnresolvedReferences
def button_select_recommended_autosomal_test():
    """ Routine to setup Recommended selections (default on entry; can be restored by user button).  """
    global checkbuttons_results, generateSelectedFilesButton, kits
    global combinedButton, g23andmev3Button, g23andmev5Button, g23andmev35Button

    # Unfortunately, no key index. So position dependent. Happens to be 0, 1, 3 and 5.
    button_deselect_every_autosomal_test()
    combinedButton.select()     ;  checkbuttons_results[kits.index("CombinedKit")] = True
    g23andmev3Button.select()   ;  checkbuttons_results[kits.index("23andMe_V3")]  = True
    g23andmev5Button.select()   ;  checkbuttons_results[kits.index("23andMe_V5")]  = True
    g23andmev35Button.select()  ;  checkbuttons_results[kits.index("23andMe_V35")] = True
    generateSelectedFilesButton.configure(state="normal")


# noinspection PyUnresolvedReferences
def button_select_every_autosomal_test():
    global checkbuttons_results, generateSelectedFilesButton
    global combinedButton
    global g23andmev3Button, g23andmev4Button, g23andmev5Button, g23andmeAPIButton, g23andmev35Button
    global ancestryv1Button, ancestryv2Button, ftdnav2Button, ftdnav3Button
    global ldnav1Button, ldnav2Button, myheritagev1Button, myheritagev2Button
    global mthfrGenButton, generaButton, meuDNAButton, g1240KButton, hov1Button, g1240hoButton

    combinedButton.select()
    g23andmev3Button.select()
    g23andmev4Button.select()
    g23andmev5Button.select()
    g23andmeAPIButton.select()
    g23andmev35Button.select()
    ancestryv1Button.select()
    ancestryv2Button.select()
    ftdnav2Button.select()
    ftdnav3Button.select()
    ldnav1Button.select()
    ldnav2Button.select()
    myheritagev1Button.select()
    myheritagev2Button.select()
    mthfrGenButton.select()
    generaButton.select()
    meuDNAButton.select()
    g1240KButton.select()
    hov1Button.select()
    g1240hoButton.select()

    for i in range(len(checkbuttons_results)):
        checkbuttons_results[i] = True
    generateSelectedFilesButton.configure(state="normal")


# noinspection PyUnresolvedReferences
def button_deselect_every_autosomal_test():
    global checkbuttons_results, generateSelectedFilesButton
    global combinedButton
    global g23andmev3Button, g23andmev4Button, g23andmev5Button, g23andmeAPIButton, g23andmev35Button
    global ancestryv1Button, ancestryv2Button, ftdnav2Button, ftdnav3Button
    global ldnav1Button, ldnav2Button, myheritagev1Button, myheritagev2Button
    global mthfrGenButton, generaButton, meuDNAButton, g1240KButton, hov1Button, g1240hoButton

    combinedButton.deselect()
    g23andmev3Button.deselect()
    g23andmev4Button.deselect()
    g23andmev5Button.deselect()
    g23andmeAPIButton.deselect()
    g23andmev35Button.deselect()
    ancestryv1Button.deselect()
    ancestryv2Button.deselect()
    ftdnav2Button.deselect()
    ftdnav3Button.deselect()
    ldnav1Button.deselect()
    ldnav2Button.deselect()
    myheritagev1Button.deselect()
    myheritagev2Button.deselect()
    mthfrGenButton.deselect()
    generaButton.deselect()
    meuDNAButton.deselect()
    g1240KButton.deselect()
    hov1Button.deselect()
    g1240hoButton.deselect()

    for i in range(len(checkbuttons_results)):
        checkbuttons_results[i] = False
    generateSelectedFilesButton.configure(state="disabled")


def adna_checkbutton_click(kit):
    global checkbuttons_results, generateSelectedFilesButton, kits

    i = kits.index(kit)
    checkbuttons_results[i] = not checkbuttons_results[i]

    number_of_selected_buttons = 0
    for i in range(len(checkbuttons_results)):
        if checkbuttons_results[i]:
            number_of_selected_buttons += 1
    generateSelectedFilesButton.configure(state="normal" if number_of_selected_buttons > 0 else "disabled")


def get_target_type_suffix(target_type_name_all):
    if "FTDNA" in target_type_name_all or "MyHeritage" in target_type_name_all:
        return ".csv"
    elif "LDNA" in target_type_name_all or "23andMe" in target_type_name_all or "Ancestry" in target_type_name_all:
        return ".txt"
    else:
        return ".err"                             # REH 20Mar2020 Refactor
    # Todo need to add MTHFR, Genera, meuDNA, 1240K and HOv1


def button_generate_selected_autosomal():
    """ Run Microarray format generator (originally based on Extract23) to generate selected formats. """
    global selectAutosomalFormatsWindow, kits

    if wgse.BAM.chrom_types['A'] == 0:
        # Need at least Autosomes present to create microarray file. Button should not have been available but ...
        # Todo should provide a pop-up warning message that doing nothing here
        cancel_autosomal_formats_window()
        return

    # Why the warning? Need to document in the manual. HG19 not good? Comparisons seem to show close enough
    # if wgse.BAM.Refgenome != "hs37d5":
    #    wgse_message("warning", 'NoHs37d5WarningTitle', True,
    #                 wgse.lang.i18n['NoHs37d5WarningText'].replace("{{refgenome}}", wgse.BAM.Refgenome))

    # selectAutosomalFormatsWindow.withdraw()

    # Call internal button to reuse or create a new CombinedKit file in the Output Diretory from the current BAM file
    reuse = _button_CombinedKit(selectAutosomalFormatsWindow)

    # By this point, we should have a .zip and .txt CombinedKit file

    out_FPB = wgse.outdir.FPB
    out_oFPB = wgse.outdir.oFPB
    out_FPB_cmbkit = f'{wgse.outdir.FPB}_CombinedKit'
    CombinedKitTXT_oFN = nativeOS(f'{out_FPB_cmbkit}.txt')
    CombinedKitZIP_oFN = nativeOS(f'{out_FPB_cmbkit}.zip')

    if os.path.exists(CombinedKitTXT_oFN) and \
        os.path.getmtime(CombinedKitTXT_oFN) > os.path.getmtime(wgse.BAM.file_oFN) and \
        os.path.getsize(CombinedKitTXT_oFN) > minCbnKitSize:
        # Now for EACH format selected (other than CombinedKit), add a run of aconv.py (our script)
        #  to take the CombinedKit content and subset it to create the respective format. Note that we
        #  are simply running aconv.py as a stand-alone python program; not importing and calling it.
        # Todo modify below to simply call aconv routine, within same python thread, and use python zip command

        python = wgse.python3x_qFN
        aconv = f'"{wgse.prog_FP}aconv.py"'
        cmdzip = wgse.zipx_qFN

        commands = ""
        stop = len(checkbuttons_results) - (0 if wgse.DEBUG_MODE else 6)    # New formats not quite ready for primetime
        for i in range(stop):
            out_compressed_oFN = f'{out_oFPB}_{kits[i]}.zip'
            if i != kits.index("CombinedKit") and checkbuttons_results[i] and \
                not (os.path.exists(out_compressed_oFN) and
                     os.path.getmtime(out_compressed_oFN) > os.path.getmtime(CombinedKitTXT_oFN) and
                     os.path.getsize(out_compressed_oFN) > minMicZipsize):
                # These files only take a minute or two each to create; but if doing all 20, takes awhile. So check
                suffix = get_target_type_suffix(kits[i])
                out_uncompressed_qFN = f'"{out_FPB}_{kits[i]}{suffix}"'
                out_compressed_qFN = f'"{out_FPB}_{kits[i]}.zip"'
                commands += (
                    f'echo "Generating microarray file for format {kits[i]}"\n'
                    f'{python} {aconv} {kits[i]} "{out_FPB_cmbkit}.txt" "{out_FPB}" "{wgse.reflib.cma_FP}"\n'
                    f'{cmdzip} -mj {out_compressed_qFN} {out_uncompressed_qFN}\n'  # Do a move; delete uncompressed
                )

        run_bash_script("ButtonMicroarrayDNA", commands, parent=selectAutosomalFormatsWindow)

    # Handle CombinedKit file cleanup based on whether requested and/or reused an existing copy
    wgse.tempf.list.append(CombinedKitTXT_oFN)      # Always delete the uncompressed version;
    # if reused or asked to save by clicking or big enough (to not be in an error), then delete the zip file also
    if not (reuse or checkbuttons_results[kits.index("CombinedKit")]) or \
       os.path.getsize(CombinedKitZIP_oFN) < minCbnKitSize:
        wgse.tempf.list.append(CombinedKitZIP_oFN)

    # button_autosomal_stats()

    cancel_autosomal_formats_window()


def _button_CombinedKit(parent_window, parallel=False):
    """
      CombinedKit merged Microarray File Creator (internal button)
      Loosely, originally based on the Extract23 script: https://github.com/tkrahn/extract23/blob/master/extract23.sh
      With NO parallelization extension shown in: https://gist.github.com/tkrahn/ef62cfaab678f447ea53ddee09ce0eb2
      Original Extract23 did this only for a single 23andMe v3 template.  Here we created a CombinedKit template
      of the merger of all known formats.  Then wrote our own aconv.py module to subset the CombinedKit for each
      desired file format.
    """

    out_FPB_cmbkit = f'{wgse.outdir.FPB}_CombinedKit'
    out_FB_cmbkit = f'{wgse.BAM.file_FB}_CombinedKit'
    out_FP = f'{wgse.outdir.FP}'

    python = wgse.python3x_qFN
    lifthg38 = f'"{wgse.prog_FP}hg38tohg19.py"'
    bamfile = f'"{wgse.BAM.file_FN}"'
    refgen_qFN = wgse.BAM.Refgenome_qFN

    # Select reference genome. Add liftover command to modify CombinedKit if GRCh38/HG38
    refVCFtab_qFN = wgse.reflib.get_reference_vcf_qFN(wgse.BAM.Build, wgse.BAM.SNTypeC, ttype="Microarray")
    if not os.path.isfile(nativeOS(unquote(refVCFtab_qFN))):
        wgse_message("error", "errRefLibMissingFileTitle", True,
                     wgse.lang.i18n["errRefLibMissingCombKitTable"].replace("{{FILE}}", refVCFtab_qFN))
        return

    liftover_tohg19 = ""        # Default; no liftover needed
    if wgse.BAM.Build == 38:
        refmod = "hg38" if wgse.BAM.SNTypeC == "Chr" else "GRCh38" if wgse.BAM.SNTypeC == "Num" else "error"
        liftover_tohg19 = f'{python} {lifthg38} "{out_FPB_cmbkit}" {refmod}\n'
        if refmod == "error":
            wgse_message("error", "errBuildUnkSeqNameTitle", True,
                         wgse.lang.i18n["errBuildUnkSeqNameType"].replace("{{SNType}}", wgse.BAM.SNTypeC))
            return
    elif wgse.BAM.Build == 99:
        wgse_message("error", "errLitfoverUnsupTitle", True,
                     wgse.lang.i18n["errLoftoverUnsuported"].replace("{{START}}", "T2T").replace("{{END}}", "19"))
        return

    # Setup additional filenames and paths
    temp_pileup_vcf = f'"{wgse.tempf.FP}CombKit_pileup.vcf.gz"'
    temp_called_vcf = f'"{wgse.tempf.FP}CombKit_called.vcf.gz"'
    temp_annotated_vcf = f'"{wgse.tempf.FP}CombKit_annotated.vcf.gz"'
    temp_result_tab = f'"{wgse.tempf.FP}CombKit_result.tab"'
    temp_sorted_result_tab = f'"{wgse.tempf.FP}CombKit_result_sorted.tab"'
    # In case DEBUG mode on and Temp directory not being cleared out; lets clear these files to make sure not reused
    wgse.tempf.list.append(temp_called_vcf) if os.path.exists(temp_called_vcf) else ""
    wgse.tempf.list.append(temp_annotated_vcf) if os.path.exists(temp_annotated_vcf) else ""
    wgse.tempf.list.append(temp_result_tab) if os.path.exists(temp_result_tab) else ""
    wgse.tempf.list.append(temp_sorted_result_tab) if os.path.exists(temp_sorted_result_tab) else ""

    temp_head_qFN = f'"{wgse.reflib.cma_FP}raw_file_templates/head/23andMe_V3.txt"'
    ploidy_qFN = f'"{wgse.reflib.cma_FP}ploidy.txt"'

    # samtools = wgse.samtoolsx_qFN
    bcftools = wgse.bcftoolsx_qFN
    tabix = wgse.tabixx_qFN
    sed = wgse.sedx_qFN
    sort = wgse.sortx_qFN
    cat = wgse.catx_qFN
    cmdzip = wgse.zipx_qFN
    cmdunzip = wgse.unzipx_qFN
    cpus = wgse.os_threads

    CombinedKitTXT_oFN = nativeOS(f'{out_FPB_cmbkit}.txt')
    CombinedKitZIP_oFN = nativeOS(f'{out_FPB_cmbkit}.zip')

    if os.path.exists(CombinedKitZIP_oFN) and \
       os.path.getmtime(CombinedKitZIP_oFN) > os.path.getmtime(wgse.BAM.file_oFN) and \
       os.path.getsize(CombinedKitZIP_oFN) > minCbnKitSize:

        if os.path.exists(CombinedKitTXT_oFN) and \
           os.path.getmtime(CombinedKitTXT_oFN) > os.path.getmtime(wgse.BAM.file_oFN) and \
           os.path.getsize(CombinedKitTXT_oFN) > minCbnKitSize:
            commands = ""
        else:  # Uncompress existing CombinedKit file in the Output directory; creating txt version
            commands = f'{cmdunzip} -oj "{out_FPB_cmbkit}.zip" -d "{out_FP}" "{out_FB_cmbkit}.txt"\n'
        reuse = True
    else:  # Have to create a CombinedKit from scratch; use Outdir for result but delete if not requested later
        if wgse.reflib.missing_refgenome(refgen_qFN, parent_window):
            return False       # missing_refgenome() reports sn error if reference genome does not exist

        '''
          Based on the (original) Extract23 Windows script: https://github.com/tkrahn/extract23/blob/master/extract23.sh
          With parallelization extension shown in: https://gist.github.com/tkrahn/ef62cfaab678f447ea53ddee09ce0eb2
          Original Extract23 did this only for a single 23andMe v3 template.  Here we created a CombinedKit template
          of the merger of all known formats.  Then wrote our own aconv.py module to subset the CombinedKit for each
          desired file format. Note that CombinedKit only includes valid calls; not no-call indication for ones missing.
        '''
        # Old versus new pileup command; old using samtools generates a warning that we should switch to new
        #   {samtools} mpileup -B    -C 50 -r {fchrom} -l {reftab_qFN} -f {refgenome} -hu {bamfile}
        #   {bcftools} mpileup -B -I -C 50 -r {fchrom} -T {reftab_qFN} -f {refgenome} -Ou {bamfile}
        # Note -B required for Nanopore long read tools; minimal effect on paired-end sequencer output
        commands = (
            f'{bcftools} mpileup -B -I -C 50 -T {refVCFtab_qFN} -f {refgen_qFN} -Ou {bamfile} -o {temp_pileup_vcf}\n'
            f'{bcftools} call {temp_pileup_vcf} --ploidy-file {ploidy_qFN} -V indels -m -P 0 --threads {cpus} -Oz -o {temp_called_vcf}\n'
        )

        # Note: liftover_hg38to19 will be Null if not a Build38 BAM so simply puts an empty line when not used
        # Todo more simplication and reduction of intermediate files possible? sed/sort/cat on one line? Avoid more
        #  compression and tabix calls by piping directly?
        commands += (
            f'{tabix} -p vcf {temp_called_vcf}\n'
            f'{bcftools} annotate -Oz -a {refVCFtab_qFN} -c CHROM,POS,ID {temp_called_vcf} > {temp_annotated_vcf}\n'
            f'{tabix} -p vcf {temp_annotated_vcf}\n'
            f'{bcftools} query -f \'%ID\t%CHROM\t%POS[\t%TGT]\n\' {temp_annotated_vcf} -o {temp_result_tab}\n'
            f'{sed} \'s/chr//; s/\tM\t/\tMT\t/g; s/\///; s/\.\.$/--/; s/TA$/AT/; s/TC$/CT/; s/TG$/GT/; s/GA$/AG/; '
            f'    s/GC$/CG/; s/CA$/AC/\' {temp_result_tab} | {sort} -t $\'\t\' -k2,3 -V > {temp_sorted_result_tab}\n'
            f'{cat} {temp_head_qFN} {temp_sorted_result_tab} > "{out_FPB_cmbkit}.txt"\n'
            f'{liftover_tohg19}'
            f'{cmdzip} -j "{out_FPB_cmbkit}.zip" "{out_FPB_cmbkit}.txt"\n'
        )

        reuse = False

    run_bash_script("ButtonCombinedKit", commands, parent=parent_window)
    # Todo really should check if completed OK; but general issue with all run_bash_script calls

    return reuse


def cancel_autosomal_formats_window():
    """ Remove Autosomal formats window in preparation for restoring WGSE Main Window """
    from mainwindow import mainwindow_resume   # Have to localize due to import loop

    global selectAutosomalFormatsWindow

    try:        # May have been destroyed before this call
        selectAutosomalFormatsWindow.destroy()
    except Exception:
        pass
    if __name__ != '__main__':      # for standalone operation
        mainwindow_resume()    # sets up for return to main window


# ***************MAIN PROGRAM (stand-alone)*******************
# Todo standalone program in development; not really been tested of late
if __name__ == '__main__':
    from mainwindow import set_BAM_file, set_output_path       # Need local import to avoid loop

    if len(sys.argv) != 3:
        print(f'Usage: microarray outdir bam_file (fully qualified paths)')
        exit()

    """ Microarray start as an independent task (main program) """
    print('Starting WGS Extract Microarray Combined Kit Generator ...')
    wgse.init(interactive=False)

    # Could potentially just get parameters from stored settings JASON file? Or make args optional and call button?
    set_output_path(sys.argv[1])
    if not (wgse.BAM or set_BAM_file(sys.argv[2])):
        exit()      # Already reported issue

    _button_CombinedKit(None)
