@echo on

docker run --rm -it -v "C:\SimLab\CrystalScatter\Mcp:/home/user/project" crystalscattercons-run sh -c "cd /home/user/project; sas_scatter2Cons --mcpinp testmcppar.txt --mcplog testdocker.log --mcpimg testdockerout.png"

pause
