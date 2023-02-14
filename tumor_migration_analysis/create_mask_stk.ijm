//macro for making a mask stack from a brightfield time series
//after running this you will be prompted to select an image
//beware that the mask image will be saved with the same name as the file, just with a "_mask" appended at the end

setBatchMode(true);

open();

dir = getInfo("image.directory");
filename = getInfo("image.filename");

if (endsWith(filename, ".tif")) {

	setOption("BlackBackground", true);
	run("Select None");
	
	//makes new image stack
	input = getImageID();
	Stack.getDimensions(width, height, no_channels, no_slices, no_frames);
	newImage("mask", "8-bit black", width, height, no_frames);
	
	//goes through image series
	for (i=1; i<no_frames+1; i++) {
		
		print("analyzing frame " + i);
		
		selectImage(input);
		setSlice(i);
		run("Duplicate...", "title=tmp");
		run("Variance...", "radius=2 stack");
		run("Auto Threshold", "method=Otsu white");
		for (j=0; j<5; j++) {
			run("Dilate", "stack");
		}
		run("Fill Holes", "stack");
		run("Connected Components Labeling", "connectivity=4 type=[16 bits]");
		selectWindow("tmp-lbl");
		run("Keep Largest Region");
		selectWindow("tmp-lbl-largest");
		run("Select All");
		run("Copy");
		selectWindow("mask");
		setSlice(i);
		run("Paste");
		run("Select None");
		selectWindow("tmp");
		close();
		selectWindow("tmp-lbl");
		close();
		selectWindow("tmp-lbl-largest");
		close();
	
	}
	
	selectWindow("mask");
	run("Grays");
	setSlice(1);
	
	saveAs("tiff",dir + "/" + replace(filename,".tif","_mask.tif"));
			
} else {
	exit("File must be a tif image");
}

close("*");
