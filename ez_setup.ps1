$mydir = Get-Location
cd ~
mkdir envs
cd envs
python -m venv indramodel
.\indramodel\Scripts\Activate.ps1
pip install wheel indra
cd $mydir
