import sys
import os
import glob
import subprocess
#Module to install required packages locally
#Written by Buys de Barbanson, 2016.
#py@buysdb.nl

def addToPath(path):

	if not path in os.environ['PATH'] and os.path.exists(path):
		os.environ['PATH'] = '%s:%s' %(path, os.environ['PATH'] )



#unfortunately pip.get_installed_distributions() misses modules added to the pypath.
#This is a routine which resolves this
def getInstalledModules():

	proc = subprocess.Popen(["pip freeze"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()

	parts = out.split('\n')
	i = []
	for part in parts:
		installedModule = part.strip().split('==')
		if len(installedModule)>=2:
			i.append(installedModule[0])
	return(i)

#Check availability of required modules, and install locally if not present.
def require( requiredModules ):
	availableModules = getInstalledModules()
	dependenciesRequired = False
	for required in requiredModules:
		if not required in availableModules:
			dependenciesRequired = True
			#print('%s is not installed on the system, trying local copy' % required)
		else:
			#print('%s OK'%required)
			pass
	#local dependency compatibility
	if dependenciesRequired:
		codeDirectory = os.path.dirname(os.path.realpath(sys.argv[0]))

		depDir = '%s/%s' %(codeDirectory, 'dependencies')
		if not os.path.exists(depDir):
			try:
				os.makedirs(depDir)
			except:
				print('Failed to create %s, please check your permissions' % depDir)
				exit()

		#Locate and include all locally installed modules:
		for eggPath in glob.glob('%s/dependencies/lib/python2.7/site-packages/*.egg' % codeDirectory):
			sys.path.insert(0,os.path.abspath(eggPath))


		#Check wether still some modules are missing:
		packageDir = depDir + '/lib/python2.7/site-packages'
		if not os.path.exists(packageDir):
			os.makedirs(packageDir)
		#Add the package directory to the pythonpath, prevents a mad easy_install.
		if not "PYTHONPATH" in os.environ:
			os.environ["PYTHONPATH"] = packageDir
		else:
			os.environ["PYTHONPATH"] = ('%s:'% packageDir )+os.environ["PYTHONPATH"]
		availableModules = getInstalledModules()

		didInstallNewPackage = False
		for required in requiredModules:
			with open('%s/installLog' % depDir, 'a') as f:
				if not required in availableModules:
					print('Running local install of %s to %s' % (required,codeDirectory))
					proc = subprocess.Popen(['easy_install --prefix %s %s' % (depDir, required)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
					(out, err) = proc.communicate()
					f.write('Installation of %s\nStandard out:\n%s\nStandard error\n%s\n'%(required, out,err))
					didInstallNewPackage = True
					#pip.main(['install','--install-option','--prefix', "%s" % depDir, required])
		#When we did install a module locally, re-configure the PATH
		if didInstallNewPackage:
			#Locate and include all locally installed modules:
			for eggPath in glob.glob('%s/dependencies/lib/python2.7/site-packages/*.egg' % codeDirectory):
				sys.path.insert(0,os.path.abspath(eggPath))

		availableModules = getInstalledModules()
		dependenciesRequired = False
		for required in requiredModules:
			if not required in availableModules:
				dependenciesRequired = True
				print('[FATAL] %s is still missing. Please check your permissions on %s'% (required,depDir))
				exit()
	print('All dependencies are present')

