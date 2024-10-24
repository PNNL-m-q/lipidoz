"""
lipidoz/gui/__main__.py

Dylan Ross (dylan.ross@pnnl.gov)

    Runs the GUI application when invoked directly
"""


# TODO: Add a couple CLI args to change hoe the app starts. In particular, 
#       it would be nice to have a --view [results.loz] option to directly 
#       load the results viewer for a specified results file without going 
#       through the startup window and file chooser. 


from lipidoz.gui.app import LOzApp


def main():

    app = LOzApp()
    app.run()


if __name__ == '__main__':
    main()
