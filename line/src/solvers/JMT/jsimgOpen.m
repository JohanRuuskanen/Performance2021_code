function ret = jsimgOpen()
cmd = ['java --illegal-access=permit -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimgraph.mainGui.JSIMGraphMain']
ret = system(cmd);
end

