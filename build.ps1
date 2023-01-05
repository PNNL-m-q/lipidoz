python -m PyInstaller --clean -w -F -y -i='lipidoz.ico' --add-data='lipidoz.ico;.' `
    --collect-submodules='sklearn' `
    --collect-data='hdf5plugin' `
    LipidOz.py
[System.Reflection.Assembly]::LoadWithPartialName('System.Windows.Forms')
[System.Windows.Forms.MessageBox]::Show('build done')
