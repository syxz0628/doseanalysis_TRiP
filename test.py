path2doselist=['c/lung.nrrd','c/lung.nrrd']
onehedfile=path2doselist[0][path2doselist[0].rfind('/')+1:]
onehedfile=onehedfile.replace('.nrrd','.hed')
print(onehedfile)