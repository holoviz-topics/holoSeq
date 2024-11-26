#!/usr/bin/env python
import argparse
import sys
from CPT_GFFParser import gffParse, gffSeqFeature



def convertSeqFeat(inFeat, defaultSource="gffSeqFeature"):
    featLoc = inFeat.location
    IDName = inFeat.id
    qualDict = inFeat.qualifiers
    parentCands = inFeat.qualifiers.get("Parent", [])
    for x in parentCands:
        if x == inFeat.id:  # Cannot allow self-loops
            raise Exception(
                "Cannot convert SeqRecord, feature %s lists itself as a parent feature"
                % (x.id)
            )
    if "codon_start" in inFeat.qualifiers.keys():
        phaseIn = int(inFeat.qualifiers["codon_start"][0])
    else:
        phaseIn = 0
    if "score" in inFeat.qualifiers.keys():
        scoreIn = float(inFeat.qualifiers["score"][0])
    else:
        scoreIn = "."
    if "source" in inFeat.qualifiers.keys():
        sourceIn = inFeat.qualifiers["source"][0]
    else:
        sourceIn = defaultSource

    return gffSeqFeature(
        featLoc,
        inFeat.type,
        "",
        featLoc.strand,
        IDName,
        qualDict,
        [],
        None,
        None,
        phaseIn,
        scoreIn,
        sourceIn,
    )


def convertSeqRec(
    inRec, defaultSource="gffSeqFeature", deriveSeqRegion=True, createMetaFeat=None
):
    # Assumes an otherwise well-constructed SeqRecord that just wants to replace its features with gffSeqFeatures
    if not isinstance(inRec, list):
        inRec = [inRec]
    outRec = []
    for rec in inRec:
        topList = []
        childList = []
        noIDList = []
        expectedParents = 0
        maxLoc = 0
        for feat in rec.features:
            if "Parent" in feat.qualifiers.keys():
                childList.append(
                    (convertSeqFeat(feat, defaultSource), [])
                )  # Possible to have more than one parent
                expectedParents += len(feat.qualifiers["Parent"])
                # lastCount += childList[-1][1]
            elif (
                feat.id and feat.id != "<unknown id>"
            ):  # Do not accept the default value
                topList.append(convertSeqFeat(feat, defaultSource))
            else:
                noIDList.append()
            maxLoc = max(maxLoc, feat.location.end)
        if deriveSeqRegion:
            rec.annotations["sequence-region"] = "%s 1 %s" % (rec.id, str(maxLoc))
        noEdit = False
        foundParCount = 0
        while not noEdit:
            noEdit = True
            for childInd in range(0, len(childList)):
                for i in childList[childInd][0].qualifiers["Parent"]:
                    nextChild = False
                    for topFeat in topList:
                        checkTree = [topFeat]
                        for parCand in checkTree:
                            nextPar = False
                            checkTree += parCand.sub_features
                            for foundPrior in childList[childInd][1]:
                                if parCand.id == foundPrior:
                                    nextPar = True
                                    break
                            if nextPar:
                                break
                            if i == parCand.id:
                                parCand.sub_features.append(childList[childInd][0])
                                childList[childInd] = (
                                    childList[childInd][0],
                                    childList[childInd][1] + [i],
                                )
                                noEdit = False
                                nextChild = True
                                foundParCount += 1
                                break
                        if nextChild:
                            break
            if noEdit and foundParCount < expectedParents:
                badFeats = ""
                for x in childList:
                    if len(x[0].qualifiers["Parent"]) != len(x[1]):
                        badFeats += x.id + ", "
                sys.stderr.write(
                    "Unable to convert SeqRecord %s: could not find parents for features [%s]\n"
                    % (rec.id, badFeats)
                )
        if createMetaFeat:
            qualDict = {}
            for x in rec.annotations.keys():
                outVal = ""
                if isinstance(rec.annotations[x], list):
                    outVal = " ".join(rec.annotations[x])
                else:
                    outVal = str(rec.annotations[x])
                outVal = outVal.replace("\n", " ")
                qualDict[x] = [outVal]
            topList.append(
                gffSeqFeature(
                    FeatureLocation(0, maxLoc),
                    createMetaFeat,
                    "",
                    0,
                    IDName,
                    qualDict,
                    [],
                    None,
                    None,
                    0,
                    ".",
                    defaultSource,
                )
            )
        topList = sorted(topList, key=lambda feature: feature.location.start)
        rec.features = topList
        outRec.append(rec)
    return outRec


def parseGFF(gff3):
    recs = gffParse(
        gff3, suppressMeta=0, pragmaPriority=True, pragmaOverridePriority=True
    )
    featList = []
    for x in recs[0].features:
        featList.append(x)
        for y in x.sub_features:
            featList.append(y)
            for z in y.sub_features:
                featList.append(z)
    print(1, len(featList))
    for x in featList:
        x.sub_features = []
    recs[0].features = featList
    recs = convertSeqRec(recs, defaultSource="ConvertSeqRec")
    featList = []
    for x in recs[0].features:
        featList.append(x)
        for y in x.sub_features:
            featList.append(y)
            for z in y.sub_features:
                featList.append(z)
    print(2, len(featList))
    gff = {"start": [], "end": [], "strand": [], "type": [], "qualifiers": []}
    for f in featList:
        gff["start"].append(int(f.location.start))
        gff["end"].append(int(f.location.end))
        gff["strand"].append(f.location.strand)
        gff["type"].append(f.type)
        gff["qualifiers"].append(f.qualifiers)

    return gff


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract specified qualifers from features in GFF3", epilog=""
    )
    parser.add_argument("gff3", type=argparse.FileType("r"), help="GFF3 File")
    args = parser.parse_args()
    gff = parseGFF(args.gff3)
    print('types=', list(set(gff['type'])))
    print('qual=',(gff['qualifiers']))

