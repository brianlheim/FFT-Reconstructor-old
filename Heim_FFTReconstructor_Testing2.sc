Heim_FFTReconstructor_Testing2 {
	classvar <>nMagSums = 32;
	classvar <>toleranceCoeff = 15;
	classvar <>sampleRate = 48e3;
	classvar <>databaseSearchIndexPercentage = 0.2;
	classvar <>analysisFunc = \totalMag;
	classvar collectedSourcePostfix = "_collected_source.wav";
	classvar reconstructionPostfix = "_recon.wav";
	classvar databaseMagsPostfix = "_dbmags";
	classvar databaseSumsPostfix = "_dbsums";
	classvar databaseOrderPostfix = "_dborder";

	var <fftsize, <hopsizeRead, <hopsizeWrite, <sourcePaths, <targetPath, <reconstructionPrefix, reconstructionPath, databaseMagsPath, databaseSumsPath, databaseOrderPath;

	var sourceFiles, collectedSourceFile, targetFile, reconstructionFile, databaseMagsFile, databaseSumsFile, databaseOrderFile;
	var costable, imagtable, window;

	var orderingArray;
	var analysisArray;
	var partialMagSumsArray;
	var cumulativeError = 0.0;
	var cumulativeTime = 0.0;

	*new {
		arg fftsize, hopsizeRead, hopsizeWrite, sourcePaths, targetPath, reconstructionPrefix;
		^super.newCopyArgs(fftsize, hopsizeRead, hopsizeWrite, sourcePaths, targetPath, reconstructionPrefix).init;
	}

	init {
		File.exists(targetPath).not.if {Error("Target file does not exist!").throw};
		costable = Signal.fftCosTable(fftsize);
		imagtable = Signal.newClear(fftsize);
		window = Signal.hanningWindow(fftsize) / (fftsize / hopsizeWrite / 2);
		targetFile = SoundFile.new(targetPath);
		if(sourcePaths.isKindOf(Array).not) {sourcePaths = sourcePaths.bubble};
		sourceFiles = sourcePaths.collect({
			|filename|
			File.exists(filename).not.if {Error("File does not exist:" + filename).throw};
			SoundFile.new(filename)
		});
		collectedSourceFile = SoundFile.new(PathName(targetPath).pathOnly ++ reconstructionPrefix ++ collectedSourcePostfix);
		reconstructionPath = PathName(targetPath).pathOnly ++ reconstructionPrefix ++ reconstructionPostfix;
		reconstructionFile = SoundFile.new(reconstructionPath);
		databaseMagsPath = PathName(targetPath).pathOnly ++ reconstructionPrefix ++ databaseMagsPostfix;
		databaseSumsPath = PathName(targetPath).pathOnly ++ reconstructionPrefix ++ databaseSumsPostfix;
		databaseOrderPath = PathName(targetPath).pathOnly ++ reconstructionPrefix ++ databaseOrderPostfix;
	}

	/////////////// database methods ////////////////

	makeDatabase {
		var flag = false, index = 0, sourceArrays = Array.fill(fftsize.div(hopsizeRead), {
			FloatArray.newClear(hopsizeRead)
		});
		var databaseMagsPath_temp = databaseMagsPath++"__temp", databaseSumsPath_temp = databaseSumsPath++"__temp";
		var databaseSumsFile_temp, databaseMagsFile_temp;
		var sumsBuffer = FloatArray.newClear(nMagSums), magsBuffer = FloatArray.newClear(fftsize.div(2)+1);

		analysisArray = FloatArray();

		"collecting source files".postln;
		this.collectSourceFiles();
		"source files collected".postln;

		"creating temp databases".post;
		databaseMagsFile_temp = File(databaseMagsPath_temp, "w");
		databaseSumsFile_temp = File(databaseSumsPath_temp, "w");
		collectedSourceFile.openRead;

		sourceArrays.do(collectedSourceFile.readData(_));
		this.makeTempDatabaseEntries(this.processData(sourceArrays.reduce('++')), databaseMagsFile_temp, databaseSumsFile_temp);

		while {flag.not} {
			((index = index + 1) % 1000 == 0).if { ".".post } {};
			sourceArrays = sourceArrays.rotate(-1);
			collectedSourceFile.readData(sourceArrays.last);

			if(sourceArrays.last.size < hopsizeRead) {
				"EOF reached.".postln;
				flag = true;
			} {
				this.makeTempDatabaseEntries(this.processData(sourceArrays.reduce('++')), databaseMagsFile_temp, databaseSumsFile_temp);
			};
		};

		collectedSourceFile.close;
		databaseMagsFile_temp.close;
		databaseSumsFile_temp.close;
		"temp databases complete".postln;

		"creating final databases".postln;
		orderingArray = Int32Array.newFrom(analysisArray.order);
		analysisArray = FloatArray.newFrom(analysisArray[orderingArray]);
		databaseMagsFile = File(databaseMagsPath, "w");
		databaseSumsFile = File(databaseSumsPath, "w");
		databaseOrderFile = File(databaseOrderPath, "w");
		databaseMagsFile_temp = File(databaseMagsPath_temp, "r");
		databaseSumsFile_temp = File(databaseSumsPath_temp, "r");

		databaseOrderFile.write(orderingArray);
		databaseOrderFile.write(analysisArray);
		orderingArray.do {
			|orderEntry, i|
			databaseMagsFile_temp.seek(fftsize.div(2)+1*4*orderEntry,0);
			databaseSumsFile_temp.seek(nMagSums*4*orderEntry,0);
			databaseMagsFile_temp.read(magsBuffer);
			databaseSumsFile_temp.read(sumsBuffer);
			databaseMagsFile.write(magsBuffer);
			databaseSumsFile.write(sumsBuffer);
		};

		databaseMagsFile.close;
		databaseOrderFile.close;
		databaseSumsFile.close;
		databaseMagsFile_temp.close;
		databaseSumsFile_temp.close;
		File.delete(databaseMagsPath_temp);
		File.delete(databaseSumsPath_temp);
		"final databases complete".postln;
	}

	collectSourceFiles {
		collectedSourceFile.headerFormat_("WAV").numChannels_(1).sampleRate_(sampleRate).sampleFormat_("int16");
		collectedSourceFile.openWrite;

		sourceFiles.do {
			|sourceFile|
			var flag = false;
			var readArray = FloatArray.newClear(hopsizeRead);

			(sourceFile.sampleRate != sampleRate).if {
				("source file has wrong sample rate:" + sourceFile.path).warn;
			};
			(sourceFile.numChannels != 1).if {
				Error("source file is not mono:" + sourceFile.path).throw;
			};
			sourceFile.openRead;

			while {flag.not} {
				sourceFile.readData(readArray);
				if(readArray.size < hopsizeRead) {
					flag = true;
					readArray = readArray.extend(hopsizeRead, 0);
				};
				collectedSourceFile.writeData(readArray);
			};

			sourceFile.close;
		};

		collectedSourceFile.close;
	}

	////////////////////////////////////////////////

	makeReconstruction {
		arg minErrorInit = 5;

		reconstructionFile.headerFormat_("WAV").numChannels_(1).sampleRate_(sampleRate).sampleFormat_("int16");
		reconstructionFile.openWrite;
		reconstructionFile = this.reconstructStartAndLoop(minErrorInit);
		reconstructionFile.close;
		"done".postln;
	}

	continueReconstruction {
		arg minErrorInit = 5;

		// create the continuation file
		var tempPathName = reconstructionPath++"__temp";
		var reconstPathName = reconstructionPath;

		// copy over everything
		var startFrame = this.copyPreviousReconstruction(reconstPathName, tempPathName);

		reconstructionFile = this.reconstructStartAndLoop(minErrorInit, startFrame - hopsizeWrite, true);
		reconstructionFile.close;
		"done".postln;
	}

	copyPreviousReconstruction {
		arg reconstPathName, tempPathName;

		var copyArray = FloatArray.newClear(hopsizeWrite);
		var tempFile, numHops;

		// copy the file at reconstruction file to the temp path
		File.copy(reconstPathName, tempPathName);
		tempFile = SoundFile.openRead(tempPathName);
		File.delete(reconstPathName);
		// create the reconstruction file and open it
		reconstructionFile = SoundFile.new(reconstPathName);
		reconstructionFile.headerFormat_("WAV").numChannels_(1).sampleRate_(sampleRate).sampleFormat_("int16");
		reconstructionFile.openWrite;
		numHops = tempFile.numFrames.div(hopsizeWrite);

		numHops.do {
			tempFile.readData(copyArray);
			reconstructionFile.writeData(copyArray);
		};

		tempFile.close;
		File.delete(tempPathName);

		^numHops * hopsizeWrite;
	}

	reconstructStartAndLoop {
		arg minErrorInit = 5, startFrame = 0, isContinuation = false;
		var targetArrays, time, targetComplex, targetMatch, constructedSig, index = 0, flag = false;
		var sigCreationArray = FloatArray.newClear(fftsize);

		targetFile.openRead;
		targetFile.seek(startFrame);
		this.createPartialMagSumsArray;
		this.createOrderArrays;
		collectedSourceFile.openRead;
		databaseMagsFile = File.new(databaseMagsPath, "r");

		targetArrays = Array.fill(fftsize.div(hopsizeWrite), {
			FloatArray.newClear(hopsizeWrite)
		});
		targetArrays.do(targetFile.readData(_));

		time = Process.elapsedTime;

		targetComplex = this.processData(targetArrays.reduce('++'));
		targetMatch = this.getDatabaseMatch(targetComplex.magnitude[0..(fftsize.div(2))], minErrorInit);
		constructedSig = this.createSig(targetMatch[0], sigCreationArray);
		minErrorInit = targetMatch[1];

		while {flag.not} {
			var tempsig;

			index.postln;
			time = this.logTime(time);

			targetArrays = targetArrays.rotate(-1);
			targetFile.readData(targetArrays.last);

			if(targetArrays.last.size < hopsizeWrite) {
				flag = true;
				// write the last bit of data
				this.writeConstructedSig(constructedSig, fftsize);
			} {
				if(isContinuation) {
					isContinuation = false; // ignore the first hop if we're doing a continuation
				} {
					this.writeConstructedSig(constructedSig, hopsizeWrite);
				};

				targetComplex = this.processData(targetArrays.reduce('++'));
				targetMatch = this.getDatabaseMatch(targetComplex.magnitude[0..(fftsize.div(2))], minErrorInit * 2);
				tempsig = this.createSig(targetMatch[0], sigCreationArray);
				minErrorInit = targetMatch[1];
				constructedSig = this.overlapSigs(constructedSig, tempsig);
			};

			index = index + 1;
		};

		targetFile.close;
		databaseMagsFile.close;
		collectedSourceFile.close;

		^reconstructionFile;
	}

	/////////////////////// database searching ////////////////////

	getDatabaseMatch {
		arg targetMags, minError;

		var match, matchIndex, finalMinError;

		// outermost loop: keep going until we found a match
		while {matchIndex.isNil} {

			match = this.searchDatabaseForMatch(targetMags, minError);
			matchIndex = match[0];
			finalMinError = match[1];

			if(matchIndex.isNil && (minError == 0)) {minError = 0.001};
			minError = minError * 2;
		};
		^[orderingArray[matchIndex], finalMinError];
	}

	searchDatabaseForMatch {
		arg targetMags, minError;

		var minErrorInit = minError, skipCount = 0, lineLength = (fftsize.div(2)+1*4);
		var targetSums = this.createPartialMagSums(targetMags), magBuffer = FloatArray.newClear(fftsize.div(2)+1);
		var tolerance = minError * toleranceCoeff; // magic
		var tempError, tempTol, matchIndex, q;
		var nearestIndex = analysisArray.indexIn(this.analysisFunc(targetMags));
		var limits = this.makeSearchIndexLimits(nearestIndex);

		format("tol: %\tmaxError: %", tolerance.round(0.0001), minError.round(0.0001)).postln;

		for(limits[0], limits[1], {
			arg magSumsIndex;
			var sourceSums = partialMagSumsArray[magSumsIndex];
			tempTol = (targetSums - sourceSums).stdDevPop(0);
			if(tempTol < tolerance) {
				databaseMagsFile.seek(lineLength * magSumsIndex, 0).read(magBuffer);
				tempError = (targetMags - magBuffer).stdDevPop(0);

				if(tempError < minError) {
					// MATCH FOUND!
					minError = tempError;
					matchIndex = magSumsIndex;
					q = tempTol / minErrorInit;
					if(minError <= 0.01) {
						this.logMatchInfo(minError, matchIndex, skipCount, q);
						^[matchIndex, minError];
					}
				};
			} {
				skipCount = skipCount + 1;
			};
		});

		if(matchIndex.notNil) {this.logMatchInfo(minError, matchIndex, skipCount, q)};

		^[matchIndex, minError];
	}

	//////////// helpers ////////////////////////

	processData {
		arg array;
		^Signal.newFrom(array).fft(imagtable, costable);
	}

	makeTempDatabaseEntries {
		arg complex, magsFile, sumsFile;
		var mag = complex.magnitude[0..fftsize.div(2)];
		var magsums = this.createPartialMagSums(mag);
		sumsFile.write(magsums);
		magsFile.write(FloatArray.newFrom(mag));
		analysisArray = analysisArray.add(this.analysisFunc(mag));
	}

	createPartialMagSums {
		arg mag;
		var size = fftsize.div(2).div(nMagSums);
		^FloatArray.fill(nMagSums, {
			|i|
			mag[(i*size+1)..((i+1)*size+1)].sum;
		});
	}

	analysisFunc {
		arg mags;
		^switch(analysisFunc)
		{\totalMag} {mags.sum}
		{\centerMagWeight} {(mags.collect(_*_).sum)/(mags.size-1)}; // ignoring DC bin
	}

	createSig {
		arg matchIndex, array;
		collectedSourceFile.seek(matchIndex * hopsizeRead, 0).readData(array);
		^Signal.newFrom(array) * window;
	}

	writeConstructedSig {
		arg signal, length;
		var data = signal[0..(length-1)];
		reconstructionFile.writeData(data);
	}

	overlapSigs {
		arg signal1, signal2;

		^if(fftsize == hopsizeWrite) {signal2} {
			for(0, fftsize-hopsizeWrite-1) {
				|i|
				signal1[i] = signal1[i+hopsizeWrite] + signal2[i];
			};
			for(fftsize-hopsizeWrite, signal1.size-1) {
				|i|
				signal1[i] = signal2[i];
			};
			signal1;
		}
	}

	logMatchInfo {
		arg minError, matchIndex, skipCount, q;
		format("minerror: %\tminindex: %", minError.round(0.0001), matchIndex).postln;
		format("magskips: %\tmatch tol / max error: %", skipCount, q.round(0.01)).postln;
		cumulativeError = cumulativeError + minError;
		format("cumulativeError: %", cumulativeError).postln;
	}

	logTime {
		arg time;

		var temptime = Process.elapsedTime;
		cumulativeTime = cumulativeTime + (temptime - time);
		format("elapsed time: %\ttotal time: %", (temptime-time).round(0.01), cumulativeTime.round(0.01)).postln;
		^temptime;
	}

	createPartialMagSumsArray {
		// read the database
		var magSumsArraySize;
		"creating partial mag sums array".post;
		databaseSumsFile = File(databaseSumsPath, "r");
		magSumsArraySize = databaseSumsFile.length.div(4).div(nMagSums);
		partialMagSumsArray = Array.fill(magSumsArraySize, {
			|i|
			var array;
			databaseSumsFile.read(array = FloatArray.newClear(nMagSums));
			if(i % 1000 == 999) {".".post};
			array;
		});

		"\nEOF reached".postln;
		databaseSumsFile.close;
	}

	createOrderArrays {
		var arraySize;
		if(orderingArray.notNil) {^nil}; // we've already created the arrays in makeDatabase

		"creating ordering array".postln;
		databaseOrderFile = File(databaseOrderPath, "r");
		arraySize = databaseOrderFile.length.div(8); // because half is in 4-byte and half is in 4-byte entries
		orderingArray = Int32Array.newClear(arraySize);
		databaseOrderFile.read(orderingArray);
		analysisArray = FloatArray.newClear(arraySize);
		databaseOrderFile.read(analysisArray);

		"EOF reached".postln;
		databaseOrderFile.close;
	}

	makeSearchIndexLimits {
		arg indexOf;
		var res = indexOf / analysisArray.size;
		res = res + (databaseSearchIndexPercentage * [-1,1]);
		res = (res * analysisArray.size).asInteger;
		res = res.clip(0,analysisArray.size-1);
		^res;
	}

}