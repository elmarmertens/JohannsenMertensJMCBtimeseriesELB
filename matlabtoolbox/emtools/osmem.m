function memused = osmem()

allvars = whos;
memused = sum([allvars.bytes]);
