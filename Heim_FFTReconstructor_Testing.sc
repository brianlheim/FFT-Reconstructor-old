Heim_FFTReconstructor_Testing {
	classvar <>nMagSums = 32;
	classvar <>toleranceCoeff = 15;
	classvar <>sampleRate = 48e3;
	classvar collectedSourcePostfix = "_collected_source.wav";
	classvar reconstructionPostfix = "_recon.wav";
	classvar databaseMagsPostfix = "_dbmags";
	classvar databaseSumsPostfix = "_dbsums";

	var <fftsize, <hopsize, <sourcePaths, <targetPath, <reconstructionPrefix, reconstructionPath, databaseMagsPath, databaseSumsPath;

	var sourceFiles, collectedSourceFile, targetFile, reconstructionFile, databaseMagsFile, databaseSumsFile;
	var costable, imagtable, window;

	var partialMagSumsArray;
	var cumulativeError = 0.0;
	var cumulativeTime = 0.0;

	*new {
		arg fftsize, hopsize, sourcePaths, targetPath, reconstructionPrefix;
		^super.newCopyArgs(fftsize, hopsize, sourcePaths, targetPath, reconstructionPrefix).init;
	}

	init {
		File.exists(targetPath).not.if {Error("Target file does not exist!").throw};
		costable = Signal.fftCosTable(fftsize);
		imagtable = Signal.newClear(fftsize);
		window = Signal.hanningWindow(fftsize) / (fftsize / hopsize / 2);
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
	}

	/////////////// database methods ////////////////

	makeDatabase {
		var flag = false, index = 0, sourceArrays = Array.fill(fftsize.div(hopsize), {
			FloatArray.newClear(hopsize)
		});

		"collecting source files".postln;
		this.collectSourceFiles();
		"source files collected".postln;

		"creating databases".post;
		databaseMagsFile = File(databaseMagsPath, "w");
		databaseSumsFile = File(databaseSumsPath, "w");
		collectedSourceFile.openRead;

		sourceArrays.do(collectedSourceFile.readData(_));
		this.writeDatabaseLines(this.processData(sourceArrays.reduce('++')));

		while {flag.not} {
			((index = index + 1) % 1000 == 0).if { ".".post } {};
			sourceArrays = sourceArrays.rotate(-1);
			collectedSourceFile.readData(sourceArrays.last);

			if(sourceArrays.last.size < hopsize) {
				"EOF reached.".postln;
				flag = true;
			} {
				this.writeDatabaseLines(this.processData(sourceArrays.reduce('++')));
			};
		};

		collectedSourceFile.close;
		databaseMagsFile.close;
		databaseSumsFile.close;
		"databases complete".postln;
	}

	collectSourceFiles {
		collectedSourceFile.headerFormat_("WAV").numChannels_(1).sampleRate_(sampleRate).sampleFormat_("int16");
		collectedSourceFile.openWrite;

		sourceFiles.do {
			|sourceFile|
			var flag = false;
			var readArray = FloatArray.newClear(hopsize);

			(sourceFile.sampleRate != sampleRate).if {
				("source file has wrong sample rate:" + sourceFile.path).warn;
			};
			(sourceFile.numChannels != 1).if {
				Error("source file is not mono:" + sourceFile.path).throw;
			};
			sourceFile.openRead;

			while {flag.not} {
				sourceFile.readData(readArray);
				if(readArray.size < hopsize) {
					flag = true;
					readArray = readArray.extend(hopsize, 0);
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
		var targetArrays, window, time, targetComplex, targetMatch, constructedSig, index, flag = false;

		// create the continuation file
		var tempPathName = reconstructionPath++"__temp";
		var reconstPathName = reconstructionPath;

		// copy over everything
		var startFrame = this.copyPreviousReconstruction(reconstPathName, tempPathName);

		reconstructionFile = this.reconstructStartAndLoop(minErrorInit, startFrame - hopsize, true);
		reconstructionFile.close;
		"done".postln;
	}

	copyPreviousReconstruction {
		arg reconstPathName, tempPathName;

		var copyArray = FloatArray.newClear(hopsize);
		var tempFile, numHops;

		// copy the file at reconstruction file to the temp path
		File.copy(reconstPathName, tempPathName);
		tempFile = SoundFile.openRead(tempPathName);
		File.delete(reconstPathName);
		// create the reconstruction file and open it
		reconstructionFile = SoundFile.new(reconstPathName);
		reconstructionFile.headerFormat_("WAV").numChannels_(1).sampleRate_(sampleRate).sampleFormat_("int16");
		reconstructionFile.openWrite;
		numHops = tempFile.numFrames.div(hopsize);

		numHops.do {
			tempFile.readData(copyArray);
			reconstructionFile.writeData(copyArray);
		};

		tempFile.close;
		File.delete(tempPathName);

		^reconstructionFile.numFrames;
	}

	reconstructStartAndLoop {
		arg minErrorInit = 5, startFrame = 0, isContinuation = false;
		var targetArrays, time, targetComplex, targetMatch, constructedSig, index = 0, flag = false;
		var sigCreationArray = FloatArray.newClear(fftsize);

		targetFile.openRead;
		targetFile.seek(startFrame);
		this.createPartialMagSumsArray;
		collectedSourceFile.openRead;
		databaseMagsFile = File.new(databaseMagsPath, "r");

		targetArrays = Array.fill(fftsize.div(hopsize), {
			FloatArray.newClear(hopsize)
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

			if(targetArrays.last.size < hopsize) {
				flag = true;
				// write the last bit of data
				this.writeConstructedSig(constructedSig, fftsize);
			} {
				if(isContinuation) {
					isContinuation = false; // ignore the first hop if we're doing a continuation
				} {
					this.writeConstructedSig(constructedSig, hopsize);
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
		^[matchIndex, finalMinError];
	}

	searchDatabaseForMatch {
		arg targetMags, minError;

		var minErrorInit = minError, skipCount = 0, lineLength = (fftsize.div(2)+1*4);
		var targetSums = this.createPartialMagSums(targetMags), magBuffer = FloatArray.newClear(fftsize.div(2)+1);
		var tolerance = minError * toleranceCoeff; // magic
		var tempError, tempTol, matchIndex, q;

		format("tol: %\tmaxError: %", tolerance.round(0.0001), minError.round(0.0001)).postln;

		partialMagSumsArray.do {
			|sourceSums, magSumsIndex|
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
		};

		if(matchIndex.notNil) {this.logMatchInfo(minError, matchIndex, skipCount, q)};

		^[matchIndex, minError];
	}

	//////////// helpers ////////////////////////

	processData {
		arg array;
		^Signal.newFrom(array).fft(imagtable, costable);
	}

	writeDatabaseLines {
		arg complex;
		var mag = complex.magnitude[0..fftsize.div(2)];
		var magsums = this.createPartialMagSums(mag);
		databaseSumsFile.write(magsums);
		databaseMagsFile.write(FloatArray.newFrom(mag));
	}

	createPartialMagSums {
		arg mag;
		var size = fftsize.div(2).div(nMagSums);
		^FloatArray.fill(nMagSums, {
			|i|
			mag[(i*size+1)..((i+1)*size+1)].sum;
		});
	}

	createSig {
		arg matchIndex, array;
		collectedSourceFile.seek(matchIndex * hopsize, 0).readData(array);
		^Signal.newFrom(array) * window;
	}

	writeConstructedSig {
		arg signal, length;
		var data = signal[0..(length-1)];
		reconstructionFile.writeData(data);
	}

	overlapSigs {
		arg signal1, signal2;

		^if(fftsize == hopsize) {signal2} {
			for(0, fftsize-hopsize-1) {
				|i|
				signal1[i] = signal1[i+hopsize] + signal2[i];
			};
			for(fftsize-hopsize, signal1.size-1) {
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


}
