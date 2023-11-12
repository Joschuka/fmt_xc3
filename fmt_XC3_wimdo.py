from inc_noesis import *
import os

global bLoadSeveralWimdos
bLoadSeveralWimdos = True
global bLoadAnims
bLoadAnims = True

def registerNoesisTypes():
    handle = noesis.register("XC3 import", ".wimdo")
    noesis.setHandlerTypeCheck(handle, CheckType)
    noesis.setHandlerLoadModel(handle, LoadModel)
    noesis.addOption(handle, "-wimdoList", "Wimdo list", noesis.OPTFLAG_WANTARG)
    noesis.addOption(handle, "-wimdoAnm", "Wimdo anim", noesis.OPTFLAG_WANTARG)
    noesis.addOption(handle, "-wimdoSkelOnly", "Wimdo skel only", 0)
    noesis.addOption(handle, "-wimdoNoAnm", "Wimdo model only", 0)
    return 1

def CheckType(data):
    bs = NoeBitStream(data)
    if bs.readUInt() == 0x4D584D44:        
        return 1
    return 0
    
def murmurHash(key): # https://github.com/flagship-io/flagship-python-sdk/blob/6f938688428cfb0e9b358dcdeb4787c5c1e719dc/flagship/helpers/murmur32x86.py
    c1 = -0x3361d2af
    c2 = 0x1b873593
    h1 = 0
    pos = 0
    end = len(key)
    k1 = 0
    # k2 = 0
    shift = 0
    # bits = 0
    nBytes = 0

    while pos < end:
        char = key[pos]
        code = ord(key[pos])
        pos += 1
        if code < 0x80:
            k2 = code
            bits = 8
        elif code < 0x800:
            k2 = (0xC0 | (code >> 6) | ((0x80 | (code & 0x3F)) << 8))
            bits = 16
        elif code < 0xD800 or code > 0xDFFF or pos >= end:
            k2 = (0xE0 | (code >> 12) | ((0x80 | (code >> 6 & 0x3F)) << 8) | ((0x80 | (code & 0x3F)) << 16))
            bits = 24
        else:
            utf32 = ord(key[pos])
            pos += 1
            utf32 = (code - 0xD7C0 << 10) + (utf32 & 0x3FF)
            k2 = (0xff & (0xF0 | (utf32 >> 18)) | ((0x80 | (utf32 >> 12 & 0x3F)) << 8) | (
                    (0x80 | (utf32 >> 6 & 0x3F)) << 16) | ((0x80 | (utf32 & 0x3F)) << 24))
            bits = 32
        k1 = k1 | (k2 << shift)
        shift += bits
        if shift >= 32:

            k1 = ((k1 * c1) & 0xFFFFFFFF) - 2 ** 32
            k1 = ((k1 << 15) & 0xFFFFFFFF) | ((k1 & 0xFFFFFFFF) >> 17)
            k1 = ((k1 * c2) & 0xFFFFFFFF) - 2 ** 32
            h1 = h1 ^ k1

            h1 = (((h1 << 13) & 0xFFFFFFFF) - 2 ** 32) | ((h1 & 0xFFFFFFFF) >> 19)
            h1 = (((h1 * 5) & 0xFFFFFFFF) - 2 ** 32) + (-0x19ab949c)
            shift -= 32
            if shift != 0:
                k1 = k2 >> (bits - shift)
            else:
                k1 = 0
            nBytes += 4

    if shift > 0:
        nBytes += (shift >> 3)

        k1 = ((k1 * c1) & 0xFFFFFFFF)
        k1 = (((k1 << 15) & 0xFFFFFFFF) - 2 ** 32) | ((k1 & 0xFFFFFFFF) >> 17)
        k1 = ((k1 * c2) & 0xFFFFFFFF)
        h1 = h1 ^ k1

    h1 = h1 ^ nBytes

    h1 = h1 ^ ((h1 & 0xFFFFFFFF) >> 16)

    h1 = ((h1 * -0x7a143595) & 0xFFFFFFFF) - 2 ** 32
    h1 = (h1 ^ ((h1 & 0xFFFFFFFF) >> 13))
    h1 = ((h1 * -0x3d4d51cb) & 0xFFFFFFFF) - 2 ** 32
    h1 = (h1 ^ ((h1 & 0xFFFFFFFF) >> 16))

    return h1 & 0xFFFFFFFF

def sample(coeffs, t):
    a,b,c,d = coeffs
    return a * (t*t*t) + b *(t*t) + c*t + d

def processV1Anim(bs, jointList, animName, nameToFinalID):
    keyFramedBoneList = []
    bs.seek(0x51)
    bIsModelSpace = bs.readByte()
    bs.seek(0x53)
    bIsAdditive = bs.readByte()
    bs.seek(0x5C)
    animLength = bs.readUInt()
    bs.seek(0x78)
    offs, entryCount, _ = bs.readUInt64(), bs.readUInt(), bs.readUInt()
    fpsInv = 0.03333334  
    
    bs.seek(offs)
    infos = []
    for i in range(entryCount):
        infos.append([[bs.readUInt64(), bs.readUInt(), bs.readUInt()] for j in range(3)])   
    
    bs.seek(0x28)
    animDataOffs = bs.readUInt()
    
    #jInfo
    bs.seek(animDataOffs + 0x20)    
    jIDOffs, jCount, pad = bs.readUInt64(), bs.readUInt(), bs.readInt()
    jNameOffs, jNameCount, pad = bs.readUInt64(), bs.readUInt(), bs.readInt()
    weirdOffs, weirdCount, pad = bs.readUInt64(), bs.readUInt(), bs.readInt()
    # assert jCount == entryCount see pc000101.mot
    
    if pad == -1:
        assert(0) #more extra stuff, not expected in type 1
    
    #joint indices
    bs.seek(jIDOffs)
    jIDToInfo = []
    for i in range(jCount):
        jIDToInfo.append(bs.readShort())
        
    #joint name for remap
    animJNames = []
    for i in range(jNameCount):
        bs.seek(jNameOffs + i*8)
        bs.seek(bs.readUInt())
        animJNames.append(bs.readString())

    for jID, infoID in enumerate(jIDToInfo):
        if infoID == -1:
            continue
        info = infos[infoID]
        if animJNames and animJNames[jID] in nameToFinalID:
            actionBone = NoeKeyFramedBone(nameToFinalID[animJNames[jID]])
        else:
            actionBone = NoeKeyFramedBone(jID)
        posNoeKeyFramedValues = []
        rotNoeKeyFramedValues = []		
        scaleNoeKeyFramedValues = []
        for i, semInfo in enumerate(info):
            semOffs, kfCount = info[i][:2]
            bs.seek(semOffs)
            kfTimes, coeffs = [], []
            
            for j in range(kfCount):
                kfTimes.append(bs.readFloat())
                for k in range(4 if i == 1 else 3):                    
                    coeffs.append([bs.readFloat() for _ in range(4)])
            
            if kfCount == 1:
                x = sample(coeffs[0], 0)
                y = sample(coeffs[1], 0)
                z = sample(coeffs[2], 0)    
                if i == 1:
                    w = sample(coeffs[3], 0)
                if i == 0:
                    posNoeKeyFramedValues.append(NoeKeyFramedValue(0, NoeVec3([x, y, z])))
                if i == 1:
                    rotNoeKeyFramedValues.append(NoeKeyFramedValue(0, NoeQuat([x, y, z,w]).transpose()))
                elif i == 2:
                    scaleNoeKeyFramedValues.append(NoeKeyFramedValue(0, NoeVec3([x, y, z])))
            else:
                kfTimes.append(animLength)
                mul = 4 if i == 1 else 3
                for j in range(kfCount):
                    t0,t1 = int(kfTimes[j]), int(kfTimes[j+1])
                    for k in range(t1-t0):
                        x = sample(coeffs[mul * j], k)
                        y = sample(coeffs[mul * j+1], k)
                        z = sample(coeffs[mul * j+2], k)    
                        if i == 1:
                            w = sample(coeffs[mul * j+3], k)
                        t = (t0 + k) * fpsInv
                        if i == 0:
                            posNoeKeyFramedValues.append(NoeKeyFramedValue(t, NoeVec3([x, y, z])))
                        if i == 1:
                            rotNoeKeyFramedValues.append(NoeKeyFramedValue(t, NoeQuat([x, y, z,w]).transpose()))
                        elif i == 2:
                            scaleNoeKeyFramedValues.append(NoeKeyFramedValue(t, NoeVec3([x, y, z])))
            
        if bIsAdditive:
            actionBone.flags |= noesis.KFBONEFLAG_ADDITIVE
        if bIsModelSpace:
            actionBone.flags |= noesis.KFBONEFLAG_MODELSPACE
        actionBone.setRotation(rotNoeKeyFramedValues, noesis.NOEKF_ROTATION_QUATERNION_4,noesis.NOEKF_INTERPOLATE_LINEAR)
        actionBone.setTranslation(posNoeKeyFramedValues, noesis.NOEKF_TRANSLATION_VECTOR_3,noesis.NOEKF_INTERPOLATE_LINEAR)
        actionBone.setScale(scaleNoeKeyFramedValues, noesis.NOEKF_SCALE_VECTOR_3,noesis.NOEKF_INTERPOLATE_LINEAR)
        
        keyFramedBoneList.append(actionBone)       
    anim = NoeKeyFramedAnim(animName, jointList, keyFramedBoneList, 30)

    return anim
    
def processV0Anim(bs, jointList, animName, nameToFinalID):
    keyFramedBoneList = []
    bs.seek(0x51)
    bIsModelSpace = bs.readByte()
    bs.seek(0x53)
    bIsAdditive = bs.readByte()
    bs.seek(0x5C)
    animLength = bs.readUInt()
    bs.seek(0x78)
    offs, entryCount, _ = bs.readUInt64(), bs.readUInt(), bs.readUInt()
    fpsInv = 0.03333334  
    
    bs.seek(0x28)
    animDataOffs = bs.readUInt()
    
    #jInfo
    bs.seek(animDataOffs + 0x20)    
    jIDOffs, jCount, pad = bs.readUInt64(), bs.readUInt(), bs.readInt()
    jNameOffs, jNameCount, pad = bs.readUInt64(), bs.readUInt(), bs.readInt()
    weirdOffs, weirdCount, pad = bs.readUInt64(), bs.readUInt(), bs.readInt()
    
    if pad == -1:
        assert(0) #more extra stuff, not expected in type 0
    
    #joint indices
    bs.seek(jIDOffs)
    jIDs = []
    for i in range(jCount):
        jIDs.append(bs.readShort())
        
    #joint name for remap
    animJNames = []
    for i in range(jNameCount):
        bs.seek(jNameOffs + i*8)
        bs.seek(bs.readUInt())
        animJNames.append(bs.readString())
    
    for k, jID in enumerate(jIDs):
        if jID == -1:
            continue
        if animJNames and animJNames[jID] in nameToFinalID:
            actionBone = NoeKeyFramedBone(nameToFinalID[animJNames[jID]])
        else:
            actionBone = NoeKeyFramedBone(jID)
        posNoeKeyFramedValues = []
        rotNoeKeyFramedValues = []		
        scaleNoeKeyFramedValues = []
        for i in range(animLength):
            bs.seek(offs + i*jCount*48 + k*48)
            pos = [bs.readFloat() for _ in range(3)]
            bs.readFloat()
            rot = [bs.readFloat() for _ in range(4)]
            scale = [bs.readFloat() for _ in range(3)]
            bs.readFloat()
            posNoeKeyFramedValues.append(NoeKeyFramedValue(i*fpsInv, NoeVec3(pos)))
            rotNoeKeyFramedValues.append(NoeKeyFramedValue(i*fpsInv, NoeQuat(rot).transpose()))
            scaleNoeKeyFramedValues.append(NoeKeyFramedValue(i*fpsInv, NoeVec3(scale)))
            
        if bIsAdditive:
            actionBone.flags |= noesis.KFBONEFLAG_ADDITIVE
        if bIsModelSpace:
            actionBone.flags |= noesis.KFBONEFLAG_MODELSPACE
        actionBone.setRotation(rotNoeKeyFramedValues, noesis.NOEKF_ROTATION_QUATERNION_4,noesis.NOEKF_INTERPOLATE_LINEAR)
        actionBone.setTranslation(posNoeKeyFramedValues, noesis.NOEKF_TRANSLATION_VECTOR_3,noesis.NOEKF_INTERPOLATE_LINEAR)
        actionBone.setScale(scaleNoeKeyFramedValues, noesis.NOEKF_SCALE_VECTOR_3,noesis.NOEKF_INTERPOLATE_LINEAR)
        
        keyFramedBoneList.append(actionBone)       
    anim = NoeKeyFramedAnim(animName, jointList, keyFramedBoneList, 30)

    return anim

def processV3Anim(bs, jointList, animName, nameToFinalID):
    keyFramedBoneList = []
    
    #These entries hold interesting info, more research needed. For now, try retrieving rootMotion tracks when relevant from these
    bs.seek(0xC)
    chunk2EntryCount = bs.readUInt()
    bs.readUInt64()
    bs.seek(bs.readUInt64())
    chunk2Entries = [bs.readUInt64() for _ in range(chunk2EntryCount)]    
    
    bs.seek(0x51)
    bIsModelSpace = bs.readByte()
    bs.seek(0x53)
    bIsAdditive = bs.readByte()
    bs.seek(0x5C)
    animLength = bs.readUInt()
    bs.seek(0x78)
    arrays = [[bs.readUInt64(), bs.readUInt(), bs.readUInt()] for j in range(4)]
    bs.seek(arrays[0][0])
    infos = [[[bs.readUInt(), bs.readUInt(), bs.readUInt()] for _ in range(3)] for k in range(arrays[0][1])]
    bs.seek(arrays[1][0])
    coeffsV = [[bs.readFloat() for _ in range(4)] for k in range(arrays[1][1])]
    bs.seek(arrays[2][0])
    coeffsR = [[bs.readFloat() for _ in range(4)] for k in range(arrays[2][1])]
    bs.seek(arrays[3][0])
    timings = [bs.readUShort() for k in range(arrays[3][1])]   
    fpsInv = 0.03333334    
    
    bs.seek(0x28)
    animDataOffs = bs.readUInt()
    
    #jInfo
    bs.seek(animDataOffs + 0x20)    
    jIDOffs = bs.readUInt64()
    jCount = bs.readUInt()
    assert jCount == arrays[0][1]
    
    #hashes
    bs.seek(animDataOffs + 0x58)
    bs.seek(bs.readUInt() + 0x20)
    hashOffset = bs.readUInt64()
    hashCount = bs.readUInt()
    assert hashCount == arrays[0][1]
    
    #root motion (mostly for cutscene anims, check D&J for obvious ex. Need more testing to see if game anims have root motion somewhere too, didn't test these much)
    bs.seek(animDataOffs + 0x50)
    rootMInfoOffs = bs.readUInt()
    bHasRootM = rootMInfoOffs > 0
    if bHasRootM:
        bs.seek(rootMInfoOffs)
        rootMArrays = [[bs.readUInt64(), bs.readUInt(), bs.readUInt()] for j in range(6)]
        posNoeKeyFramedValues = []
        rotNoeKeyFramedValues = []		
        scaleNoeKeyFramedValues = []
        rootPosVals = []
        rootRotVals = []
        rootScaleVals = []
        rootTimings = []
        actionBone = NoeKeyFramedBone(0)
        
        #pos
        bs.seek(rootMArrays[0][0])
        for i in range(rootMArrays[0][1]):
            rootPosVals.append(NoeVec3([bs.readFloat(), bs.readFloat(), bs.readFloat()]))
            bs.readFloat()
        #rot
        bs.seek(rootMArrays[1][0])
        for i in range(rootMArrays[1][1]):
            rootRotVals.append(NoeQuat([bs.readFloat(), bs.readFloat(), bs.readFloat(), bs.readFloat()]).transpose())
        #scale
        bs.seek(rootMArrays[2][0])
        for i in range(rootMArrays[2][1]):
            rootScaleVals.append(NoeVec3([bs.readFloat(), bs.readFloat(), bs.readFloat()]))
            bs.readFloat()
         
        #pos kf
        bs.seek(rootMArrays[3][0])
        for k in range(rootMArrays[3][1]):
            id = bs.readUShort()
            posNoeKeyFramedValues.append(NoeKeyFramedValue(k* fpsInv, rootPosVals[id]))
            
        #rot kf
        bs.seek(rootMArrays[4][0])
        for k in range(rootMArrays[4][1]):
            id = bs.readUShort()
            rotNoeKeyFramedValues.append(NoeKeyFramedValue(k* fpsInv, rootRotVals[id]))
            
        #scale kf
        bs.seek(rootMArrays[5][0])
        for k in range(rootMArrays[5][1]):
            id = bs.readUShort()
            scaleNoeKeyFramedValues.append(NoeKeyFramedValue(k* fpsInv, rootScaleVals[id]))
            
        if bIsAdditive:
            actionBone.flags |= noesis.KFBONEFLAG_ADDITIVE
        if bIsModelSpace:
            actionBone.flags |= noesis.KFBONEFLAG_MODELSPACE
        actionBone.setRotation(rotNoeKeyFramedValues, noesis.NOEKF_ROTATION_QUATERNION_4,noesis.NOEKF_INTERPOLATE_LINEAR)
        actionBone.setTranslation(posNoeKeyFramedValues, noesis.NOEKF_TRANSLATION_VECTOR_3,noesis.NOEKF_INTERPOLATE_LINEAR)
        actionBone.setScale(scaleNoeKeyFramedValues, noesis.NOEKF_SCALE_VECTOR_3,noesis.NOEKF_INTERPOLATE_LINEAR)
        keyFramedBoneList.append(actionBone)
    else:#No "easy" rootM trakcs, try retrieving it from chunk2 entries, if it exists
        for e in chunk2Entries:
            bs.seek(e)
            candidateOffs = bs.readUInt64()
            if bs.readUInt() == animLength and bs.readUInt64() == 0xFFFFFFFFFFFFFFFF and bs.readUInt() == 0xFFFFFFFF: #figure out how it truly works eventually. Rot tracks too somewhere?
                bs.seek(candidateOffs)
                rootPosVals = [bs.read('<4f')[:-1] for i in range(animLength)]
                posNoeKeyFramedValues = [NoeKeyFramedValue(k* fpsInv, rootPosVals[k]) for k in range(animLength)]
                actionBone = NoeKeyFramedBone(0)
                actionBone.setTranslation(posNoeKeyFramedValues, noesis.NOEKF_TRANSLATION_VECTOR_3,noesis.NOEKF_INTERPOLATE_LINEAR)
                keyFramedBoneList.append(actionBone)
                bHasRootM = True
    
    #joint indices
    bs.seek(jIDOffs)
    jIDs = []
    for _ in range(jCount):
        jIDs.append(bs.readUShort())
    
    #hashes
    bs.seek(hashOffset)
    hashMap = {}
    remap = {}
    for joint in jointList:
        hashMap[murmurHash(joint.name)] = [joint.name, joint.index]
    for jID in jIDs:
        hash = bs.readUInt()
        if hash not in hashMap:
            continue #Investigate this later instead of continuing, weapon related joints?
        remap[jID] = hashMap[hash][1]
    
    for jID, info in zip(jIDs, infos):
        if jID == 0 and bHasRootM:
            continue 
        if jID not in remap:
            continue #same as above
        jID = remap[jID]
        posInfo, rotInfo, scaleInfo = info
        posNoeKeyFramedValues = []
        rotNoeKeyFramedValues = []		
        scaleNoeKeyFramedValues = []
        actionBone = NoeKeyFramedBone(jID)
        for i, semInfo in enumerate(info):
            s, trackID, e = semInfo
            kfCount = e-s
            if kfCount==1:
                if i == 1:
                    x = sample(coeffsR[trackID], 0)
                    y = sample(coeffsR[trackID + 1], 0)
                    z = sample(coeffsR[trackID + 2], 0)                
                    w = sample(coeffsR[trackID + 3], 0)
                else:
                    x = sample(coeffsV[trackID], 0)
                    y = sample(coeffsV[trackID + 1], 0)
                    z = sample(coeffsV[trackID + 2], 0)  
                if i == 0:
                    posNoeKeyFramedValues.append(NoeKeyFramedValue(0, NoeVec3([x, y, z])))
                if i == 1:
                    rotNoeKeyFramedValues.append(NoeKeyFramedValue(0, NoeQuat([x, y, z,w]).transpose()))
                elif i == 2:
                    scaleNoeKeyFramedValues.append(NoeKeyFramedValue(0, NoeVec3([x, y, z])))
            else:
                kfTimes = timings[s:e]
                kfTimes.append(animLength)
                for j in range(kfCount):
                    t0,t1 = kfTimes[j], kfTimes[j+1]
                    for k in range(t1-t0):
                        if i == 1:
                            x = sample(coeffsR[trackID + 4 * j], k)
                            y = sample(coeffsR[trackID + 4 * j + 1], k)
                            z = sample(coeffsR[trackID + 4 * j + 2], k)                        
                            w = sample(coeffsR[trackID + 4 * j + 3], k)
                        else:
                            x = sample(coeffsV[trackID + 3 * j], k)
                            y = sample(coeffsV[trackID + 3 * j + 1], k)
                            z = sample(coeffsV[trackID + 3 * j + 2], k) 
                        t = (t0 + k) * fpsInv
                        if i == 0:
                            posNoeKeyFramedValues.append(NoeKeyFramedValue(t, NoeVec3([x, y, z])))
                        if i == 1:
                            rotNoeKeyFramedValues.append(NoeKeyFramedValue(t, NoeQuat([x, y, z,w]).transpose()))
                        elif i == 2:
                            scaleNoeKeyFramedValues.append(NoeKeyFramedValue(t, NoeVec3([x, y, z])))
        if bIsAdditive:
            actionBone.flags |= noesis.KFBONEFLAG_ADDITIVE
        actionBone.setRotation(rotNoeKeyFramedValues, noesis.NOEKF_ROTATION_QUATERNION_4,noesis.NOEKF_INTERPOLATE_LINEAR)
        actionBone.setTranslation(posNoeKeyFramedValues, noesis.NOEKF_TRANSLATION_VECTOR_3,noesis.NOEKF_INTERPOLATE_LINEAR)
        actionBone.setScale(scaleNoeKeyFramedValues, noesis.NOEKF_SCALE_VECTOR_3,noesis.NOEKF_INTERPOLATE_LINEAR)
        
        keyFramedBoneList.append(actionBone)       
    anim = NoeKeyFramedAnim(animName, jointList, keyFramedBoneList, 30)

    return anim
    
def processAnims(bs, jointList, animName, nameToFinalID):
    bs.seek(0x50)
    opcode = bs.readUByte()
    if opcode == 1:
        return processV1Anim(bs, jointList, animName, nameToFinalID)
    elif opcode == 3:
       return processV3Anim(bs, jointList, animName, nameToFinalID)
    elif opcode == 0:
        return processV0Anim(bs, jointList, animName, nameToFinalID)
    
    
    

def LoadAnms(jointList, animationList, motPaths, nameToFinalID):  
    for path in motPaths:
        motData = rapi.loadIntoByteArray(path)
        if path.endswith(".mot"):
            bs2 = NoeBitStream(motData)
            bDEHack = False
            if bs2.readUInt() == 0x31636278: #compressed mot, see DE
                bs2.seek(8)
                decompSize = bs2.readUInt()
                compSize = bs2.readUInt()
                bs2.readBytes(0x20)
                motData = rapi.decompInflate(bs2.readBytes(compSize),decompSize)
                bs2 = NoeBitStream(motData)
                bDEHack = True
            bs2.seek(0xC)
            n = bs2.readUInt()
            offs = bs2.readUInt()
            bs2.seek(offs)
            offsets, sizes, names = [], [], []
            for i in range(n):
                next = bs2.tell() + 0x40
                offsets.append(bs2.readUInt())
                sizes.append(bs2.readUInt())
                bs2.readUInt()
                names.append(bs2.readString())
                bs2.seek(next)
            for offs, size, name in zip(offsets, sizes, names):
                if name.endswith(".anm") or bDEHack:
                    an = processAnims(NoeBitStream(motData[offs:offs+size]), jointList, name, nameToFinalID)
                    if an is not None:
                        animationList.append(an)
        else:
            bs2 = NoeBitStream(motData)
            lName = rapi.getLocalFileName(path)
            an = processAnims(bs2, jointList, lName, nameToFinalID)
            if an is not None:
                animationList.append(an)      
                
def LoadSkel(bs, offs):
    jointList = []
    
    bs.seek(offs + bs.readUInt() + 0x20)
    arrays = [[bs.readUInt64(), bs.readUInt(), bs.readUInt()] for j in range(8)]
    
    #parenting
    pList = []
    bs.seek(arrays[0][0] + offs)
    for i in range(arrays[0][1]):
        pList.append(bs.readShort())
    
    #names
    jNames = []
    for i in range(arrays[1][1]):
        bs.seek(arrays[1][0] + offs + i*0x10)        
        bs.seek(bs.readUInt64() + offs)
        jNames.append(bs.readString())
        
    #mats
    jMats = []
    bs.seek(arrays[2][0] + offs)
    for i in range(arrays[2][1]):
        pos = NoeVec3.fromBytes(bs.readBytes(0xC))
        bs.readUInt()
        quat = NoeQuat.fromBytes(bs.readBytes(0x10)).transpose()
        bs.read('4f')
        mat = quat.toMat43()
        mat[3] = pos
        jMats.append(mat)
        
    for i, (p, n, m) in enumerate(zip(pList, jNames, jMats)):
        jointList.append(NoeBone(i,n,m,None,p))
    jointList = rapi.multiplyBones(jointList)
    
    return jointList

def findAndLoadChr(inputName, seenChr):
    chrName = rapi.getLocalFileName(inputName)
    filesL = []
    for root, dirs, files in os.walk(rapi.getDirForFilePath(inputName)):
        for fileName in files:
            filesL.append([root, fileName])
    originalName = chrName + ".wimdo"
    discardedCandidates = {}
    while len(chrName)>=4:
        for f in filesL:                
            if f[1].startswith(chrName) and (f[1].endswith(".chr") or f[1].endswith(".arc")) and f[1] <= originalName and f[1] not in discardedCandidates:
                #did we already load that one for a previous wimdo
                if f[1] in seenChr:
                    return []
                # check if the candidate does actually have a skel:
                bs = NoeBitStream((rapi.loadIntoByteArray(os.path.join(f[0], f[1]))))
                bs.seek(0xC)
                entryCount, entryOffs = bs.readUInt(), bs.readUInt()
                bs.seek(entryOffs)
                offsets, sizes, names = [], [], []
                for i in range(entryCount):
                    next = bs.tell() + 0x40
                    offsets.append(bs.readUInt())
                    sizes.append(bs.readUInt())
                    bs.readUInt()
                    names.append(bs.readString())
                    bs.seek(next)
                for offs in offsets:
                    bs.seek(offs + 0x20)
                    if bs.readUInt64() == 5495881740129927174: #6 SKEL
                        seenChr[f[1]] = True
                        print("Loaded: " + f[1])
                        return LoadSkel(bs, offs)
                discardedCandidates[f[1]] = True
                   
        chrName = chrName[:-1]
    return None

def LoadModel(data, mdlList):
    
    global bLoadSeveralWimdos
    global bLoadAnims
    
    ctx = rapi.rpgCreateContext()
    
    wimdoPaths = [rapi.getInputName()]
    motPaths = []
    
    if noesis.optWasInvoked("-wimdoList"):
        bLoadSeveralWimdos = False            
        for p in noesis.optGetArg("-wimdoList").split(","):
            if p.endswith(".wimdo"):
                wimdoPaths.append(os.path.dirname(rapi.getInputName()) + os.sep + p)    

    if noesis.optWasInvoked("-wimdoAnm"):
        bLoadAnims = False
        for p in noesis.optGetArg("-wimdoAnm").split("||"):
            motPaths.append(p)
    
    if noesis.optWasInvoked("-wimdoNoAnm"):
        bLoadAnims = False
    
    #should probably cache paths only and check for integrity later. Oh well
    if bLoadSeveralWimdos:
        wimdoLoop = True
        while wimdoLoop:
            wimdoPath = rapi.loadPairedFileGetPath("wimdo file", ".wimdo")
            if wimdoPath is not None:
                wimdoPaths.append(wimdoPath[1])
            else:
                wimdoLoop = False
                
    
    if bLoadAnims:
        motLoop = True
        while motLoop:
            wimdoPath = rapi.loadPairedFileGetPath("Motion data", ".mot;*.anm")
            if wimdoPath is not None:
                motPaths.append(wimdoPath[1])
            else:
                motLoop = False
                
    jointList = []
    nameToFinalID = {}
    nameToFinalID[""] = -1
    seenChr = {}
    
    #Skel (chr)
    for wimdoID, wimdoPath in enumerate(wimdoPaths):
        inputName = rapi.getExtensionlessName(wimdoPath)
        jList = findAndLoadChr(inputName, seenChr)
        if jList is not None and jList:
            jList = list(jList)
            for joint in jList:
                if joint.name not in nameToFinalID:
                    joint.parentName = jList[joint.parentIndex].name if joint.parentIndex >= 0 else ""
                    if joint.parentName not in nameToFinalID:
                        assert 0
                    # Monolith typo  ch01021011
                    if joint.name == "FC_tootnbase":
                        joint.name = "FC_toothbase"
                    nameToFinalID[joint.name] = len(jointList)
                    jointList.append(NoeBone(len(jointList), joint.name, joint.getMatrix(), joint.parentName))
        elif jList is None:
            assert 0
     
    for wimdoID, wimdoPath in enumerate(wimdoPaths): 
        if noesis.optWasInvoked("-wimdoSkelOnly"):
            continue
        # Buffers (wismt)
        print(wimdoPath)
        inputName = rapi.getExtensionlessName(wimdoPath)
        if not rapi.checkFileExists(inputName+".wismt"):
            return 0
            
        bs = NoeBitStream((rapi.loadIntoByteArray(inputName+".wismt")))
        bs.seek(0xC)
        resOffs = bs.readUInt()
        bs.seek(resOffs + 0x8)
        streamEntryCount = bs.readUInt()
        streamEntryOffs = bs.readUInt()
        internalStreamCount = bs.readUInt()
        internalStreamOffs = bs.readUInt()
        
        bs.seek(streamEntryOffs + resOffs)
        streamEntries = []
        for i in range(streamEntryCount):
            streamEntries.append([bs.readUInt(), bs.readUInt(), bs.readUShort(), bs.readUShort()]) # offs, size, streamIdx, entryType
            bs.readUInt64()
        internalStreamInfo = []    
        bs.seek(internalStreamOffs + resOffs)
        for i in range(internalStreamCount):
            internalStreamInfo.append([bs.readUInt(), bs.readUInt(), bs.readUInt()])
        
        modelStream = None
        for sEntry in streamEntries:
            if sEntry[3] == 0: #model
                bs.seek(internalStreamInfo[sEntry[2]][2])
                if sEntry[0] != sEntry[1]:
                    bs.seek(8, 1)
                    decompSize = bs.readUInt()
                    compSize = bs.readUInt()
                    bs.readBytes(0x20)
                    modelStream = rapi.decompInflate(bs.readBytes(compSize),decompSize)
                    
        bs = NoeBitStream(modelStream)
        vBuffOfs, vBuffCount, idxBuffOfs, idxBuffCount = bs.readUInt(), bs.readUInt(), bs.readUInt(), bs.readUInt()
        bs.read('6i')
        morphOffs, buffSize, buffOffs, _, skinBufferInfoOfs = bs.readUInt(), bs.readUInt(), bs.readUInt(), bs.readUInt(), bs.readUInt()
        
        # vBuffs 
        bs.seek(vBuffOfs)
        vBufferInfo = [[bs.readUInt() for i in range(8)] for _ in range(vBuffCount)]
        for vBuff in vBufferInfo:
            bs.seek(vBuff[3])
            vBuff.append([[bs.readUShort(), bs.readUShort()] for _ in range(vBuff[4])])
        # idxBuffs 
        bs.seek(idxBuffOfs)
        idxBufferInfo = [[bs.readUInt() for i in range(5)] for _ in range(idxBuffCount)]  
        wPals = [] #investigate this stuff properly at some point, see Ghidra's bookmarks for the involved functions and skinflags.Use 
        # skinBuffer
        if skinBufferInfoOfs:
            bs.seek(skinBufferInfoOfs)
            skinBufferInfo = [bs.readUInt(), bs.readUInt(), bs.readUShort(), bs.readUShort(), bs.readUInt(), bs.readUInt()]
            bs.seek(skinBufferInfo[1])
            for i in range(skinBufferInfo[0]):
                wPals.append([bs.readUInt() for _ in range(7)] + [bs.readUByte() for _ in range(4)] + [bs.readUInt(), bs.readUInt()])
        # morphs
        buffIdxToMorphBuffer = {} 

        if morphOffs:
            bs.seek(morphOffs)
            morphCount = bs.readUInt()
            morphSubOffs = bs.readUInt()
            bs.readUInt()
            morphSubOffs2 = bs.readUInt()
            for m in range(morphCount):
                bs.seek(morphSubOffs + m*0x14)
                vBufferIndex = bs.readUInt()
                offs = bs.readUInt()
                bs.seek(morphSubOffs2 + 0x10 * offs)
                defaultOfs = buffOffs+bs.readUInt()
                defaultVCount = bs.readUInt()
                defaultStride = bs.readUInt()                
                bs.seek(defaultOfs)
                defaultBuffer = bs.readBytes(defaultVCount*defaultStride)
                buffIdxToMorphBuffer[vBufferIndex] = [defaultBuffer,defaultStride]  
        
        # wimdo
        bs = NoeBitStream(data if not wimdoID else rapi.loadIntoByteArray(wimdoPath))
        bs.seek(0x8)
        modelOffs = bs.readUInt()
        # Model info
        bs.seek(modelOffs)
        bs.readUInt()
        bs.read('6f')
        meshOffs = bs.readUInt()
        meshCount = bs.readUInt()
        bs.readUInt()
        skinOffs = bs.readUInt() + modelOffs
        bs.read('21f')
        morphOffs = bs.readUInt()
        
        # Mesh
        primList = []
        for m in range(meshCount):        
            bs.seek(modelOffs + meshOffs + 0x44*m)
            primOffs = bs.readUInt()
            primCount = bs.readUInt()
            bs.read('8f')
            # Prim
            bs.seek(modelOffs + primOffs)
            for p in range(primCount):            
                renderFlags = bs.readUInt()
                skinFlags = bs.readUInt() #See Ghidra's bookmarked functions. Field used to set some globabl and local flags used in a two level lookup to grab the final w palette.
                vBufferID = bs.readUShort()
                idxBufferID = bs.readUShort()
                bs.read('4i')
                bs.readUShort()
                lod = bs.readUShort()
                bs.read('4i')
                primList.append([vBufferID, idxBufferID, lod, renderFlags, skinFlags])
        
        # skinInfo
        doJointList = []
        localIDToName = {}
        
        bs.seek(skinOffs + 0x4)
        doJointCount = bs.readUInt()
        doNameOffs = bs.readUInt()
        doJointMatOffs = bs.readUInt()
        bs.read('3i')
        doNonASJointIndicesOffs = bs.readUInt() #Not sure
        ASSectionOffs1 = bs.readUInt()
        bs.readUInt()
        ASSectionOffs = bs.readUInt()
        
        doJointNames = []
        for i in range(doJointCount):
            bs.seek(doNameOffs + skinOffs + i*0x18)
            bs.seek(skinOffs + bs.readUInt())
            doJointNames.append(bs.readString())
            localIDToName[i] = doJointNames[-1]
            
        doJointMats = []
        bs.seek(skinOffs + doJointMatOffs)
        for i in range(doJointCount):
            doJointMats.append(NoeMat44.fromBytes(bs.readBytes(0x40)).toMat43().inverse())
        
        if doNonASJointIndicesOffs: #Hack?
            if ASSectionOffs1:
                bs.seek(skinOffs + ASSectionOffs1)
                ASDataOffs, ASDataCount = bs.readUInt(), bs.readUInt()
                bs.seek(skinOffs + ASDataOffs)
                for i in range(ASDataCount):
                    bs.readUInt()
                    jID, pID = bs.readShort(), bs.readShort()
                    doJointList.append(NoeBone(jID,doJointNames[jID],doJointMats[jID],doJointNames[pID],-1))
                    bs.readBytes(0x1C)
                doJointList.sort(key=lambda x: x.index)
            
            if ASSectionOffs:
                bs.seek(skinOffs + ASSectionOffs)
                ASDataOffs, ASDataCount = bs.readUInt(), bs.readUInt()
                bs.seek(skinOffs + ASDataOffs)
                for i in range(ASDataCount):
                    jID, pID = bs.readShort(), bs.readShort()
                    doJointList.append(NoeBone(jID,doJointNames[jID],doJointMats[jID],doJointNames[pID],-1))
                    bs.readBytes(0x4C)        
                doJointList.sort(key=lambda x: x.index)
        for doJ in doJointList:
            doJ.index = len(jointList)
            if doJ.name not in nameToFinalID:
                jointList.append(doJ)
                nameToFinalID[jointList[-1].name] = jointList[-1].index
                
        jMap = []
        for u in range(doJointCount):
            doJName = localIDToName[u]
            if doJName not in nameToFinalID:
                print(doJName + " not found in map")
                doJ = NoeBone(len(jointList),doJointNames[u],doJointMats[u],None,0)
                jointList.append(doJ)
                nameToFinalID[doJName] = jointList[-1].index
            jMap.append(nameToFinalID[doJName])
        
        bs = NoeBitStream(modelStream)
        seenIdxBuffIDs = {}
        ofs, size, stride = vBufferInfo[skinBufferInfo[2]][:3]
        bs.seek(buffOffs + ofs)
        sBuf = bs.readBytes(size * stride)
        for k, prim in enumerate(primList):
            vBufID, idxBufID, lod, drawFlags, skinFlags = prim
            if lod != 1 and lod != 0:
                continue
            bRender = False
            
            ofs, size, stride = vBufferInfo[vBufID][:3]
            bs.seek(buffOffs + ofs)
            vBuf = bs.readBytes(size * stride)
            semOfs = 0
            for desc in vBufferInfo[vBufID][-1]:
                # position
                if vBufID in buffIdxToMorphBuffer:
                    mBuf, mStride = buffIdxToMorphBuffer[vBufID]
                    rapi.rpgBindPositionBufferOfs(mBuf, noesis.RPGEODATA_FLOAT,  mStride, 0)
                    rapi.rpgBindNormalBufferOfs(mBuf, noesis.RPGEODATA_UBYTE, mStride, 0xC)     
                    bRender = True
                if desc[0] == 0:
                    bRender = True                
                    rapi.rpgBindPositionBufferOfs(vBuf, noesis.RPGEODATA_FLOAT, stride, semOfs)
                # normals
                elif desc[0] == 28:
                    rapi.rpgBindNormalBufferOfs(vBuf, noesis.RPGEODATA_BYTE, stride, semOfs)
                # UV1
                elif desc[0] == 5:
                    rapi.rpgBindUV1BufferOfs(vBuf, noesis.RPGEODATA_FLOAT, stride, semOfs)
                # UV2
                elif desc[0] == 6:
                    rapi.rpgBindUV2BufferOfs(vBuf, noesis.RPGEODATA_FLOAT, stride, semOfs)
                # UV3
                elif desc[0] == 7:
                    rapi.rpgBindUVXBufferOfs(vBuf, noesis.RPGEODATA_FLOAT, stride, 2,2,semOfs)
                
                # Do skinning eventually properly, good enough for most XC3 for now.
                elif desc[0] == 3: 
                    bs2 = NoeBitStream(vBuf)
                    bs2.seek(semOfs)
                    finalSkinBuff = bytes()
                    for _ in range(size):
                        id = bs2.readUInt()
                        finalSkinBuff += sBuf[id*12:id*12 + 12]
                        if _ != size -1:
                            bs2.seek(stride-4,1)
                    rapi.rpgBindBoneIndexBufferOfs(finalSkinBuff,noesis.RPGEODATA_UBYTE,12,8,4)
                    rapi.rpgBindBoneWeightBufferOfs(finalSkinBuff,noesis.RPGEODATA_USHORT,12,0,4)
                    rapi.rpgSetBoneMap(jMap)
              
                semOfs += desc[1]
            if bRender and not idxBufID in seenIdxBuffIDs:
                rapi.rpgSetName("mesh_" + str(k))
                seenIdxBuffIDs[idxBufID] = True
                idxBOffs, idxCount = idxBufferInfo[idxBufID][:2]            
                bs.seek(buffOffs + idxBOffs)            
                rapi.rpgCommitTriangles(bs.readBytes(idxCount * 2), noesis.RPGEODATA_USHORT,idxCount, noesis.RPGEO_TRIANGLE, 1)
            rapi.rpgClearBufferBinds()
            
    animationList = []
    if jointList and motPaths:
        LoadAnms(jointList, animationList, motPaths, nameToFinalID)    
    try:
        mdl = rapi.rpgConstructModel()
    except:
        mdl = NoeModel()
    if jointList:        
        mdl.setBones(jointList)
        if animationList:
            mdl.setAnims(animationList)
            if not bLoadSeveralWimdos:
                rapi.setPreviewOption("setSkelToShow", str(1))
    mdlList.append(mdl)
    return 1