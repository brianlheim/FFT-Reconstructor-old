Heim_FFTReconstructor {
	classvar <>nMagSums = 32;
	classvar <>toleranceCoeff = 15;

	var <fftsize, <hopsize, <sourcePaths, <targetPath, <reconstructionPath, <databasePath;

	var sourceFiles, targetFile, reconstructionFile, <databaseFile;
	var costable, imagtable, window;

	var partialMagSumsArray;
	var cumulativeError = 0.0;
	var cumulativeTime = 0.0;

	*new {
		arg fftsize, hopsize, sourcePaths, targetPath, reconstructionPath, databasePath;
		^super.newCopyArgs(fftsize, hopsize, sourcePaths, targetPath, reconstructionPath, databasePath).init;
	}

	init {
		costable = Signal.fftCosTable(fftsize);
		imagtable = Signal.newClear(fftsize);
		window = Signal.hanningWindow(fftsize) / (fftsize / hopsize / 2);
		targetFile = SoundFile.new(targetPath);
		reconstructionFile = SoundFile.new(reconstructionPath);
		if(sourcePaths.isKindOf(Array).not) {sourcePaths = sourcePaths.bubble};
		sourceFiles = sourcePaths.collect(SoundFile.new(_));
	}

	processData {
		arg array;
		^Signal.newFrom(array).fft(imagtable, costable);
	}

	createPartialMagSums {
		arg mag;
		var size = fftsize.div(2).div(nMagSums);
		^FloatArray.fill(nMagSums, {
			|i|
			mag[(i*size+1)..((i+1)*size+1)].sum;
		});
	}

	writeDatabaseLine {
		arg complex;
		var mag = complex.magnitude[0..fftsize.div(2)];
		var phase = complex.phase[0..fftsize.div(2)];
		var magsums = this.createPartialMagSums(mag);
		databaseFile.write(magsums);
		databaseFile.write(FloatArray.newFrom(mag));
		databaseFile.write(FloatArray.newFrom(phase));
	}

	getDatabaseMatch {
		arg targetMags, minError;

		var match, matchMags, matchPhases, finalMinError;

		// outermost loop: keep going until we found a match
		while {matchMags.isNil} {

			match = this.searchDatabaseForMatch(targetMags, minError);
			matchMags = match[0];
			matchPhases = match[1];
			finalMinError = match[2];

			if(matchMags.isNil && (minError == 0)) {minError = 0.001};
			minError = (minError * 2);
		};
		^[matchMags, matchPhases, finalMinError];
	}

	searchDatabaseForMatch {
		arg targetMags, minError;

		var minErrorInit = minError;
		var eofFlag = false;
		var targetSums = this.createPartialMagSums(targetMags);
		var tempError, tolerance = minError * toleranceCoeff; // magic
		var matchIndex, matchMags, matchPhases;
		var skipCountMag = 0, lineLength = (fftsize.div(2)+1*8) + (nMagSums*4);
		var index = 0;
		var tempTol, q = nil;
		var magBuffer = FloatArray.newClear(fftsize.div(2)+1), phaseBuffer = FloatArray.newClear(fftsize.div(2)+1);

		format("tol: %\tmaxError: %", tolerance.round(0.0001), minError.round(0.0001)).postln;

		partialMagSumsArray.do {
			|sourceSums, magSumsIndex|
			tempTol = (targetSums - sourceSums).stdDevPop(0);
			if(eofFlag.not && (tempTol < tolerance)) {
				databaseFile.seek((lineLength * magSumsIndex) + (nMagSums*4), 0);
				databaseFile.read(magBuffer);
				tempError = (targetMags - magBuffer).stdDevPop(0);

				if(tempError < minError) {
					// MATCH FOUND!
					databaseFile.read(phaseBuffer);
					matchMags = magBuffer.copy;
					matchPhases = phaseBuffer.copy;
					minError = tempError;
					matchIndex = index;
					q = tempTol / minErrorInit;
					if(minError <= 0.01) {

						format("minerror: %\tminindex: %", minError.round(0.0001), matchIndex).postln;
						format("magskips: %\tmatch tol / max error: %", skipCountMag, q.round(0.01)).postln;
						cumulativeError = cumulativeError + minError;
						format("cumulativeError: %", cumulativeError).postln;
						^[matchMags, matchPhases, minError];
					}
				};
			} {
				skipCountMag = skipCountMag + 1;
			};
			index = index + 1;
		};
		if(matchMags.isNil.not) {
			format("minerror: %\tminindex: %", minError.round(0.0001), matchIndex).postln;
			format("magskips: %\tmatch tol / max error: %", skipCountMag, q.round(0.01)).postln;
			cumulativeError = cumulativeError + minError;
			format("cumulativeError: %", cumulativeError).postln;
		} {};
		^[matchMags, matchPhases, minError];
	}

	createSig {
		arg matchMags, matchPhases;
		var polar, imag, real;

		polar = Polar(matchMags,matchPhases);
		real = Signal.newFrom(polar.real.mirror1);
		imag = Signal.newFrom(polar.imag++(polar.imag.drop(-1).reverse.drop(-1).neg));

		^real.ifft(imag,costable).real * window;
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

	writeConstructedSig {
		arg file, signal, length;

		var data = signal[0..(length-1)];
		file.writeData(data);
	}

	makeDatabase {
		databaseFile = File(databasePath, "w");

		sourceFiles.do {
			|sourceFile|
			var eofFlag = false, index;
			var sourceArrays = Array.fill(fftsize.div(hopsize), {
				FloatArray.newClear(hopsize)
			});

			format("now reading %", sourceFile.path).post;
			sourceFile.openRead;
			sourceArrays.do(sourceFile.readData(_));

			this.writeDatabaseLine(this.processData(sourceArrays.reduce('++')));

			index = 0;
			while {eofFlag.not} {
				(index % 1000 == 0).if { ".".post } {};
				sourceArrays = sourceArrays.rotate(-1);
				sourceFile.readData(sourceArrays.last);

				if(sourceArrays.last.size < hopsize) {
					"EOF reached.".postln;
					eofFlag = true;
				} {
					this.writeDatabaseLine(this.processData(sourceArrays.reduce('++')));
				};

				index = index + 1;
			};

			sourceFile.close;
		};

		"database complete".postln;
		databaseFile.close;
	}

	reconstructStartAndLoop {
		arg theFile, minErrorInit = 5, startFrame = 0, isContinuation = false;
		var targetArrays, time, targetComplex, targetMatch, constructedSig, index, eofFlag = false;

		targetFile.openRead;
		targetFile.seek(startFrame);
		this.createPartialMagSumsArray;
		databaseFile = File.new(databasePath, "r");

		targetArrays = Array.fill(fftsize.div(hopsize), {
			FloatArray.newClear(hopsize)
		});
		targetArrays.do(targetFile.readData(_));

		time = Process.elapsedTime;

		targetComplex = this.processData(targetArrays.reduce('++'));
		targetMatch = this.getDatabaseMatch(targetComplex.magnitude[0..(fftsize.div(2))], minErrorInit);
		constructedSig = this.createSig(targetMatch[0], targetMatch[1], window);
		minErrorInit = targetMatch[2] * 2;

		index = 0;
		while {eofFlag.not} {
			var tempsig, temptime;

			index.postln;

			temptime = Process.elapsedTime;
			cumulativeTime = cumulativeTime + (temptime - time);
			format("elapsed time: %\ttotal time: %", (temptime-time).round(0.01), cumulativeTime.round(0.01)).postln;
			time = temptime;

			targetArrays = targetArrays.rotate(-1);
			targetFile.readData(targetArrays.last);

			if(targetArrays.last.size < hopsize) {
				eofFlag = true;
				// write the last bit of data
				this.writeConstructedSig(theFile, constructedSig, fftsize);
			} {
				if(isContinuation) {
					isContinuation = false; // ignore the first hop if we're doing a continuation
				} {
					this.writeConstructedSig(theFile, constructedSig, hopsize);
				};

				targetComplex = this.processData(targetArrays.reduce('++'));
				targetMatch = this.getDatabaseMatch(targetComplex.magnitude[0..(fftsize.div(2))], minErrorInit);
				tempsig = this.createSig(targetMatch[0], targetMatch[1], window);
				minErrorInit = targetMatch[2] * 2;
				constructedSig = this.overlapSigs(constructedSig, tempsig);
			};

			index = index + 1;
		};

		targetFile.close;
		databaseFile.close;

		^theFile;
	}

	createPartialMagSumsArray {
		// read the database
		var flag = false, index = 0, seekLen = (fftsize.div(2)+1*4*2);
		partialMagSumsArray = Array.new;
		databaseFile = File(databasePath, "r");

		"creating partial mag sums array".post;

		while {flag.not} {
			var buffer = FloatArray.newClear(nMagSums);
			databaseFile.read(buffer);
			if(buffer.size > 0) {
				partialMagSumsArray = partialMagSumsArray.add(buffer);
			} {
				flag = true;
			};
			databaseFile.seek(seekLen, 1);
			index = index + 1;
			if(index % 1000 == 0) {".".post};
		};
		"\nOEF reached".postln;
		databaseFile.close;
	}



	makeReconstruction {
		arg minErrorInit = 5;

		reconstructionFile.headerFormat_("WAV");
		reconstructionFile.numChannels_(1);
		reconstructionFile.sampleRate_(48e3);
		reconstructionFile.sampleFormat("int16");
		reconstructionFile.openWrite;

		reconstructionFile = this.reconstructStartAndLoop(reconstructionFile, minErrorInit);

		reconstructionFile.close;
		"done".postln;
	}

	copyPreviousReconstruction {
		arg continuePath;

		var continueFile = SoundFile(continuePath);
		var copyArray = FloatArray.newClear(hopsize);
		var numHops;

		reconstructionFile.openRead;
		numHops = reconstructionFile.numFrames.div(hopsize);
		continueFile.headerFormat_(reconstructionFile.headerFormat);
		continueFile.numChannels_(reconstructionFile.numChannels);
		continueFile.sampleRate_(reconstructionFile.sampleRate);
		continueFile.sampleFormat_(reconstructionFile.sampleFormat);
		continueFile.openWrite;

		numHops.do {
			reconstructionFile.readData(copyArray);
			continueFile.writeData(copyArray);
		};

		reconstructionFile.close;

		^[continueFile, reconstructionFile.numFrames];
	}


	continueReconstruction {
		arg minErrorInit = 5;
		var targetArrays, window, time, targetComplex, targetMatch, constructedSig, index, eofFlag = false;

		// create the continuation file
		var reconstPath = PathName(reconstructionPath);
		var continuePath = reconstPath.pathOnly ++ reconstPath.fileNameWithoutExtension ++ "_cont." ++ reconstPath.extension;

		// copy over everything
		var continuation = this.copyPreviousReconstruction(continuePath);
		var continueFile = continuation[0];
		var startFrame = continuation[1];

		continueFile.seek(startFrame);

		continueFile = this.reconstructStartAndLoop(continueFile, minErrorInit, startFrame - hopsize, true);

		continueFile.close;
		"done".postln;
	}

}
