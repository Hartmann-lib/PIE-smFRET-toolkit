import numpy as np
import os

def read_data(fileIN, boolINFO):

    # Modified Read_PTU.py script to read HT3 formats of HydraHarp and TimeHarp devices
    # This is demo code. Use at your own risk. No warranties.
    # Keno Goertz, PicoQUant GmbH, February 2018
    # Modified by Andreas Hartmann, B CUBE, TU Dresden, July 2025

    import time
    import struct
    import io

    # Tag Types
    tyEmpty8      = struct.unpack(">i", bytes.fromhex("FFFF0008"))[0]
    tyBool8       = struct.unpack(">i", bytes.fromhex("00000008"))[0]
    tyInt8        = struct.unpack(">i", bytes.fromhex("10000008"))[0]
    tyBitSet64    = struct.unpack(">i", bytes.fromhex("11000008"))[0]
    tyColor8      = struct.unpack(">i", bytes.fromhex("12000008"))[0]
    tyFloat8      = struct.unpack(">i", bytes.fromhex("20000008"))[0]
    tyTDateTime   = struct.unpack(">i", bytes.fromhex("21000008"))[0]
    tyFloat8Array = struct.unpack(">i", bytes.fromhex("2001FFFF"))[0]
    tyAnsiString  = struct.unpack(">i", bytes.fromhex("4001FFFF"))[0]
    tyWideString  = struct.unpack(">i", bytes.fromhex("4002FFFF"))[0]
    tyBinaryBlob  = struct.unpack(">i", bytes.fromhex("FFFFFFFF"))[0]

    # Record types
    rtPicoHarp300T3  = struct.unpack(">i", bytes.fromhex('00010303'))[0]
    rtPicoHarp300T2  = struct.unpack(">i", bytes.fromhex('00010203'))[0]
    rtHydraHarpT3    = struct.unpack(">i", bytes.fromhex('00010304'))[0]
    rtHydraHarpT2    = struct.unpack(">i", bytes.fromhex('00010204'))[0]
    rtHydraHarp2T3   = struct.unpack(">i", bytes.fromhex('01010304'))[0]
    rtHydraHarp2T2   = struct.unpack(">i", bytes.fromhex('01010204'))[0]
    rtTimeHarp260NT3 = struct.unpack(">i", bytes.fromhex('00010305'))[0]
    rtTimeHarp260NT2 = struct.unpack(">i", bytes.fromhex('00010205'))[0]
    rtTimeHarp260PT3 = struct.unpack(">i", bytes.fromhex('00010306'))[0]
    rtTimeHarp260PT2 = struct.unpack(">i", bytes.fromhex('00010206'))[0]
    rtGenericT3      = struct.unpack(">i", bytes.fromhex('00010307'))[0] # MultiHarpXXX and PicoHarp330
    rtGenericT2      = struct.unpack(">i", bytes.fromhex('00010207'))[0] # MultiHarpXXX and PicoHarp330

    # global variables
    global inputfile
    global outputfile
    global recNum
    global oflcorrection
    global truensync
    global dlen
    global isT2
    global globRes
    global numRecords
    global unit

    pathname = os.path.dirname(fileIN)
    filename_woEXT = os.path.splitext(os.path.basename(fileIN))[0]

    inputfile = open(fileIN, "rb")

    if boolINFO==1:
        fileOUT = os.path.join(pathname, filename_woEXT + "_INFO.txt")
        outputfile = io.open(fileOUT, "w+", encoding="utf-8")
    else:
        outputfile=''

    # Check if inputfile is a valid PTU file
    # Python strings don't have terminating NULL characters, so they're stripped
    magic = inputfile.read(8).decode("utf-8").strip('\0')
    if magic != "PQTTTR":
        print("ERROR: Magic invalid, this is not a PTU file.")
        inputfile.close()
        if boolINFO == 1: outputfile.close()
        exit(0)

    version = inputfile.read(8).decode("utf-8").strip('\0')

    if boolINFO == 1: outputfile.write("Tag version: %s\n" % version)

    # Write the header data to outputfile and also save it in memory.
    # There's no do ... while in Python, so an if statement inside the while loop
    # breaks out of it
    tagDataList = []    # Contains tuples of (tagName, tagValue)
    while True:
        tagIdent = inputfile.read(32).decode("utf-8").strip('\0')
        tagIdx = struct.unpack("<i", inputfile.read(4))[0]
        tagTyp = struct.unpack("<i", inputfile.read(4))[0]
        if tagIdx > -1:
            evalName = tagIdent + '(' + str(tagIdx) + ')'
        else:
            evalName = tagIdent
        if boolINFO == 1: outputfile.write("\n%-40s" % evalName)
        if tagTyp == tyEmpty8:
            inputfile.read(8)
            if boolINFO == 1: outputfile.write("<empty Tag>")
            tagDataList.append((evalName, "<empty Tag>"))
        elif tagTyp == tyBool8:
            tagInt = struct.unpack("<q", inputfile.read(8))[0]
            if tagInt == 0:
                if boolINFO == 1: outputfile.write("False")
                tagDataList.append((evalName, "False"))
            else:
                if boolINFO == 1: outputfile.write("True")
                tagDataList.append((evalName, "True"))
        elif tagTyp == tyInt8:
            tagInt = struct.unpack("<q", inputfile.read(8))[0]
            if boolINFO == 1: outputfile.write("%d" % tagInt)
            tagDataList.append((evalName, tagInt))
        elif tagTyp == tyBitSet64:
            tagInt = struct.unpack("<q", inputfile.read(8))[0]
            if boolINFO == 1: outputfile.write("{0:#0{1}x}".format(tagInt,18))
            tagDataList.append((evalName, tagInt))
        elif tagTyp == tyColor8:
            tagInt = struct.unpack("<q", inputfile.read(8))[0]
            if boolINFO == 1: outputfile.write("{0:#0{1}x}".format(tagInt,18))
            tagDataList.append((evalName, tagInt))
        elif tagTyp == tyFloat8:
            tagFloat = struct.unpack("<d", inputfile.read(8))[0]
            if boolINFO == 1: outputfile.write("%-3E" % tagFloat)
            tagDataList.append((evalName, tagFloat))
        elif tagTyp == tyFloat8Array:
            tagInt = struct.unpack("<q", inputfile.read(8))[0]
            if boolINFO == 1: outputfile.write("<Float array with %d entries>" % tagInt/8)
            tagDataList.append((evalName, tagInt))
        elif tagTyp == tyTDateTime:
            tagFloat = struct.unpack("<d", inputfile.read(8))[0]
            tagTime = int((tagFloat - 25569) * 86400)
            tagTime = time.gmtime(tagTime)
            if boolINFO == 1: outputfile.write(time.strftime("%a %b %d %H:%M:%S %Y", tagTime))
            tagDataList.append((evalName, tagTime))
        elif tagTyp == tyAnsiString:
            tagInt = struct.unpack("<q", inputfile.read(8))[0]
            tmp_bytes = inputfile.read(tagInt)
            try:
                tagString = tmp_bytes.decode('utf-8').strip("\0")
            except UnicodeDecodeError:
                tagString = tmp_bytes.decode('latin1','ignore').strip("\0")
            if boolINFO == 1: outputfile.write("%s" % tagString)
            tagDataList.append((evalName, tagString))
        elif tagTyp == tyWideString:
            tagInt = struct.unpack("<q", inputfile.read(8))[0]
            tagString = inputfile.read(tagInt).decode("utf-16le", errors="ignore").strip("\0")
            if boolINFO == 1: outputfile.write(tagString)
            tagDataList.append((evalName, tagString))
        elif tagTyp == tyBinaryBlob:
            tagInt = struct.unpack("<q", inputfile.read(8))[0]
            inputfile.seek(tagInt, 1)
            if boolINFO == 1: outputfile.write("<Binary blob with %d bytes>" % tagInt)
            tagDataList.append((evalName, tagInt))
        else:
            print("ERROR: Unknown tag type")
            exit(0)
        if tagIdent == "Header_End":
            break

    # Reformat the saved data for easier access
    tagNames = [tagDataList[i][0] for i in range(0, len(tagDataList))]
    tagValues = [tagDataList[i][1] for i in range(0, len(tagDataList))]

    # get important variables from headers
    numRecords = tagValues[tagNames.index("TTResult_NumberOfRecords")]
    globRes = tagValues[tagNames.index("MeasDesc_GlobalResolution")]
    binRes = tagValues[tagNames.index("MeasDesc_Resolution")]

    def readHT3(version):

        global inputfile, numRecords

        T3WRAPAROUND = 1024

        try:
            T3Record = np.fromfile(inputfile, dtype=np.uint32, count=numRecords)
        except:
            print("Problems reading the file format.")
            exit(0)

        nsync = T3Record & 0x3FF  # bitand of T3Record with 1023
        dtime = (T3Record >> 10) & 0x7FFF  # bitshift of T3Record to the right by 10 bits and bitand with 32767
        channel = (T3Record >> 25) & 0x3F  # bitshift of T3Record to the right by 25 bits and bitand with 63
        special = (T3Record >> 31) & 0x1  # bitshift of T3Record to the right by 31 bits and bitand with 1

        oldOFL = (special != 0) & (channel == 63) & (nsync == 0) # old overflow

        OFL = (special != 0) & (channel == 63) & (nsync != 0) # new overflow

        numOFL = OFL.astype(np.uint32) * nsync + oldOFL.astype(np.uint32) # number of overflows for every entry

        calcOFL = numOFL * T3WRAPAROUND # overflow delay for every entry

        OFLCORR = np.cumsum(calcOFL) # cumulative integration of overflow delay

        # truetime for every entry
        true_nSync = OFLCORR + nsync
        truetime = true_nSync * globRes * 1e9

        # final data
        mask = special == 0

        RawData = np.column_stack([channel[mask] + 1, dtime[mask], truetime[mask]])

        return RawData

    oflcorrection = 0
    dlen = 0
    if boolINFO == 1: outputfile.write("\n-----------------------\n")
    recordType = tagValues[tagNames.index("TTResultFormat_TTTRRecType")]
    if recordType == rtHydraHarpT3:
        isT2 = False
        print("HydraHarp V1 T3 data")
        unit = 'ns'
        rawData = readHT3(1)
    elif recordType == rtHydraHarp2T3:
        isT2 = False
        print("HydraHarp V2 T3 data")
        unit = 'ns'
        rawData = readHT3(2)
    elif recordType == rtTimeHarp260NT3:
        isT2 = False
        print("TimeHarp260N T3 data")
        unit = 'ns'
        rawData = readHT3(2)
    elif recordType == rtTimeHarp260PT3:
        isT2 = False
        print("TimeHarp260P T3 data")
        unit = 'ns'
        rawData = readHT3(2)
    elif recordType == rtGenericT3:
        isT2 = False
        print("PQ Generic T3 data")
        unit = 'ns'
        rawData = readHT3(2)
    else:
        print("ERROR: This script reads only T3 data of the following devices: HydraHarp, HydraHarp2, TimeHarp260N, TimeHarp260P, and generic T3")
        exit(0)

    print("%d records ..." % numRecords)

    inputfile.close()
    if boolINFO == 1: outputfile.close()

    return rawData, unit, globRes, binRes
