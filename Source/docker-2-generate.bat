@echo on

rem *** prepare the runner image
docker run -t --name scattercons-tmp  -v "%cd%:/home/user/project" g76r/qt6-runner:qt-6.8.3-debug sh -c "apt install -y libxkbcommon-dev libgl-dev; cp /home/user/project/build/sas_scatter2Cons /bin"

rem *** remove the 'old' destination image
docker image rm crystalscattercons-run

rem *** generate the destination image
docker container commit -m "ConsExec" scattercons-tmp crystalscattercons-run

rem *** remove the runner image (not used anymore)
docker container rm scattercons-tmp

pause
