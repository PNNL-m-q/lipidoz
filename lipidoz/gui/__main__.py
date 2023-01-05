"""
lipidoz/gui/__main__.py

Dylan Ross (dylan.ross@pnnl.gov)

    Runs the GUI application when invoked directly
"""


from lipidoz.gui.app import LOzApp


def main():

    app = LOzApp()
    app.run()


if __name__ == '__main__':
    main()
