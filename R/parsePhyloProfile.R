#' Get list of pre-installed NCBI taxon names
#' @description Get all NCBI taxon names from
#' "PhyloProfile/data/taxonNamesReduced.txt"
#' @export
#' @return List of taxon IDs, their full names, taxonomy ranks and parent IDs
#' obtained from "PhyloProfile/data/taxonNamesReduced.txt"
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

getNameListCr <- function() {
    nameReducedFile <- paste(
        system.file(package="PhyloCellulase"),
        "PhyloProfile/data/taxonNamesReduced.txt",
        sep="/"
    )

    if (!file.exists(nameReducedFile)) {
        utils::data(taxonNamesReduced)
    } else {
        taxonNamesReduced <- utils::read.table(
            nameReducedFile, sep = "\t", header = TRUE, fill = TRUE,
            comment.char = ""
        )
    }

    taxonNamesReduced$fullName <- as.character(taxonNamesReduced$fullName)
    taxonNamesReduced$rank <- as.character(taxonNamesReduced$rank)
    taxonNamesReduced <- taxonNamesReduced[!duplicated(taxonNamesReduced), ]

    return(taxonNamesReduced)
}

#' Get taxonomy matrix
#' @description Get the (full or subset) taxonomy matrix from
#' "data/taxonomyMatrix.txt" based on an input taxon list
#' @export
#' @param subsetTaxaCheck TRUE/FALSE subset taxonomy matrix based on input taxon
#' IDs. Default = FALSE.
#' @param taxonIDs list of input taxon IDs (e.g. ncbi1234). Default = NULL.
#' @return Data frame contains the (subset of) taxonomy matrix for list of
#' input taxa.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

getTaxonomyMatrixCr <- function(subsetTaxaCheck = FALSE, taxonIDs = NULL){
    taxonomyMatrixFile <- paste(
        system.file(package="PhyloCellulase"),
        "PhyloProfile/data/taxonomyMatrix.txt",
        sep="/"
    )

    if (!file.exists(taxonomyMatrixFile)) {
        utils::data(taxonomyMatrix)
    } else {
        taxonomyMatrix <- utils::read.table(
            taxonomyMatrixFile, sep = "\t", header = TRUE,
            stringsAsFactors = TRUE
        )
    }

    if (subsetTaxaCheck) {
        if (missing(taxonIDs)) return(taxonomyMatrix)
        taxonomyMatrix <- taxonomyMatrix[
            taxonomyMatrix$abbrName  %in% taxonIDs, ]
    }
    return(taxonomyMatrix)
}

#' Get NCBI taxon names for a selected list of taxa
#' @description Get NCBI taxon names from
#' "PhyloProfile/data/taxonNamesReduced.txt" for a list of input taxa
#' @param rankName taxonomy rank (e.g. "species","phylum",...)
#' @param taxonIDs list of taxon IDs (e.g. ncbi1234). Default = NULL
#' @return Data frame contains a list of full names, taxonomy ranks and parent
#' IDs for the input taxa.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export

getInputTaxaNameCr <- function(rankName, taxonIDs = NULL){
    # check input parameters
    if (missing(rankName)) return("No taxonomy rank name given!")
    allMainRanks <- PhyloProfile::getTaxonomyRanks()
    if (!(rankName[1] %in% allMainRanks)) return("Invalid taxonomy rank given!")
    # load list of unsorted taxa
    Dt <- getTaxonomyMatrixCr(TRUE, taxonIDs)
    # load list of taxon name
    nameList <- getNameListCr()
    # return
    choice <- data.frame(
        "ncbiID" = unlist(Dt[rankName]), stringsAsFactors = FALSE
    )
    choice <- merge(choice, nameList, by = "ncbiID", all = FALSE)
    return(choice)
}

#' Get a subset of input taxa based on a selected taxonomy rank
#' @description Get a subset of taxon ncbi IDs and names from an input list of
#' taxa based on a selected supertaxon (identified by its taxonomy rank and
#' supertaxon name or supertaxon ID).
#' @usage getSelectedTaxonNamesCr(inputTaxonIDs, rank, higherRank, higherID,
#'     higherName)
#' @param inputTaxonIDs list of input taxon IDs (e.g. c("10116", "122586"))
#' @param rank taxonomy rank of input taxa (e.g. "species")
#' @param higherRank selected taxonomy rank (e.g. "phylum")
#' @param higherID supertaxon ID (e.g. 7711). NOTE: either supertaxon ID or
#' name is required, not neccessary to give both.
#' @param higherName supertaxon name (e.g. "Chordata"). NOTE: either
#' supertaxon ID or name is required, not neccessary to give both.
#' @export
#' @return A data frame contains ncbi IDs and names of taxa from the input taxon
#' list that belong to the selected supertaxon.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

getSelectedTaxonNamesCr <- function(
    inputTaxonIDs = NULL, rank = NULL,
    higherRank = NULL, higherID = NULL, higherName = NULL
) {
    rankName <- NULL
    if (is.null(inputTaxonIDs) | is.null(rank))
        stop("Input taxa and taxonomy rank cannot be NULL!")
    taxDf <- getTaxonomyMatrixCr(TRUE, paste0("ncbi", inputTaxonIDs))
    if (is.null(higherID) & is.null(higherName))
        return(
            data.frame(
                ncbiID = taxDf$ncbiID[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                name = taxDf$fullName[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                stringsAsFactors = FALSE))
    if (is.null(higherRank)) {
        return(
            data.frame(
                ncbiID = taxDf$ncbiID[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                name = taxDf$fullName[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                stringsAsFactors = FALSE))
    } else {
        if (!is.null(higherName) & is.null(higherID)) {
            taxaList <- getNameListCr()
            superID <- taxaList$ncbiID[
                taxaList$fullName == higherName
                & taxaList$rank %in% c(higherRank, "norank")]
            customizedtaxaID <- levels(
                as.factor(taxDf[rank][taxDf[higherRank] == superID, ]))
            return(
                data.frame(
                    ncbiID = taxaList$ncbiID[
                        taxaList$rank %in% c(rank, "norank")
                        & taxaList$ncbiID %in% customizedtaxaID],
                    name = taxaList$fullName[
                        taxaList$rank %in% c(rank, "norank")
                        & taxaList$ncbiID %in% customizedtaxaID],
                    stringsAsFactors = FALSE))
        } else if (!is.null(higherID)) {
            return(
                data.frame(
                    ncbiID = taxDf$ncbiID[taxDf[,higherRank] == higherID],
                    name = taxDf$fullName[taxDf[,higherRank] == higherID],
                    stringsAsFactors = FALSE))
        }
    }
}

#' taxa2dist
#' @param x taxa matrix
#' @param varstep var-step
#' @param check check
#' @param labels labels
#' @return a distance matrix
#' @author function from taxize library

taxa2dist <- function(x, varstep = FALSE, check = TRUE, labels) {
    rich <- apply(x, 2, function(taxa) length(unique(taxa)))
    S <- nrow(x)
    if (check) {
        keep <- rich < S & rich > 1
        rich <- rich[keep]
        x <- x[, keep]
    }
    i <- rev(order(rich))
    x <- x[, i]
    rich <- rich[i]
    if (varstep) {
        add <- -diff(c(nrow(x), rich, 1))
        add <- add/c(S, rich)
        add <- add/sum(add) * 100
    }
    else {
        add <- rep(100/(ncol(x) + check), ncol(x) + check)
    }
    if (!is.null(names(add)))
        names(add) <- c("Base", names(add)[-length(add)])
    if (!check)
        add <- c(0, add)
    out <- matrix(add[1], nrow(x), nrow(x))
    for (i in seq_len(ncol(x))) {
        out <- out + add[i + 1] * outer(x[, i], x[, i], "!=")
    }
    out <- stats::as.dist(out)
    attr(out, "method") <- "taxa2dist"
    attr(out, "steps") <- add
    if (missing(labels)) {
        attr(out, "Labels") <- rownames(x)
    }
    else {
        if (length(labels) != nrow(x))
            warning("Labels are wrong: needed ", nrow(x), " got ",
                    length(labels))
        attr(out, "Labels") <- as.character(labels)
    }
    if (!check && any(out <= 0))
        warning("you used 'check=FALSE' and some distances are zero
                -- was this intended?")
    out
}

#' Create rooted tree from a taxonomy matrix
#' @export
#' @param df data frame contains taxonomy matrix used for generating tree
#' (see distDf in example)
#' @param rootTaxon taxon used for rooting the taxonomy tree
#' @importFrom ape as.phylo
#' @importFrom ape root
#' @return A rooted taxonomy tree as an object of class "phylo".
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

createRootedTreeCr <- function(df, rootTaxon = NULL){
    if (missing(df)) return("No taxonomy matrix given!")
    # calculate distance matrix
    taxdis <- tryCatch(taxa2dist(df), error = function(e) e)
    # create tree
    tree <- ape::as.phylo(stats::hclust(taxdis))
    # root tree
    if (missing(rootTaxon)) rootTaxon = tree$tip.label[1]
    if (!(rootTaxon %in% tree$tip.label)) rootTaxon = tree$tip.label[1]
    tree <- ape::root(tree, outgroup = rootTaxon, resolve.root = TRUE)
    # return
    return(tree)
}


#' Sort list of (super)taxa based on a selected reference (super)taxon
#' @usage sortInputTaxaCr(taxonIDs = NULL, rankName, refTaxon = NULL,
#'     taxaTree = NULL)
#' @param taxonIDs list of taxon IDs (e.g.: ncbi1234, ncbi9999, ...). Default =
#' NULL.
#' @param rankName working taxonomy rank (e.g. "species", "phylum",...)
#' @param refTaxon selected reference taxon. Default = NULL.
#' @param taxaTree taxonomy tree for the input taxa (optional). Default = NULL.
#' @return A taxonomy matrix for the input taxa ordered by the selected
#' reference taxon. This matrix is sorted either based on the NCBI taxonomy
#' info, or based on an user-defined taxonomy tree (if provided).
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export

sortInputTaxaCr <- function(
    taxonIDs = NULL, rankName, refTaxon = NULL, taxaTree = NULL
){
    ncbiID <- fullName <- abbrName <- NULL
    if (missing(rankName)) return("No taxonomy rank name given!")
    allMainRanks <- PhyloProfile::getTaxonomyRanks()
    if (!(rankName[1] %in% allMainRanks)) return("Invalid taxonomy rank given!")
    if (is.null(refTaxon))  refTaxon <- taxonNames$fullName[1]
    # get list of taxon names
    fullnameList <- getNameListCr()
    taxonNames <- getInputTaxaNameCr(rankName, taxonIDs)
    # get selected supertaxon ID(s)
    rankNameTMP <- taxonNames$rank[taxonNames$fullName == refTaxon]
    if (rankName == "strain") {
        superID <- fullnameList$ncbiID[fullnameList$fullName == refTaxon]
    } else
        superID <- fullnameList$ncbiID[
            fullnameList$fullName == refTaxon
            & fullnameList$rank == rankNameTMP[1]]
    # get full taxonomy data & representative taxon
    Dt <- getTaxonomyMatrixCr()
    repTaxon <- Dt[Dt[, rankName] == superID, ][1, ]
    # THEN, SORT TAXON LIST BASED ON TAXONOMY TREE
    if (is.null(taxaTree)) {
        distDf <- subset(Dt, select = -c(ncbiID, fullName))
        row.names(distDf) <- distDf$abbrName
        distDf <- distDf[, -1]
        taxaTree <- createRootedTreeCr(
            distDf, as.character(repTaxon$abbrName)
        )
    } else
        taxaTree <- ape::root(
            taxaTree,outgroup=as.character(repTaxon$abbrName),resolve.root=TRUE)
    taxonList <- PhyloProfile::sortTaxaFromTree(taxaTree)
    sortedDt <- Dt[match(taxonList, Dt$abbrName), ]
    # subset to get list of input taxa only
    sortedDt <- subset(sortedDt, abbrName %in% taxonIDs)
    # get only taxonIDs list of selected rank and rename columns
    sortedOut <- subset(
        sortedDt,
        select = c("abbrName", "ncbiID", "fullName", as.character(rankName)))
    colnames(sortedOut) <- c("abbrName", "species", "fullName", "ncbiID")
    # add name of supertaxa into sortedOut list
    sortedOut <- merge(
        sortedOut, fullnameList, by = "ncbiID", all.x = TRUE, sort = FALSE)
    sortedOut$species <- paste0("ncbi", sortedOut$species)
    ## create new column for sorted supertaxon
    indexSpec <- unlist(lapply(
        seq_len(nlevels(as.factor(sortedOut$fullName.y))),function (x) 1000+x))
    indexSpecDf <- data.frame(
        fullName.y = unique(as.character(sortedOut$fullName.y)),
        sortedSupertaxon = paste0(
            indexSpec, "_", unique(as.character(sortedOut$fullName.y))
        ), stringsAsFactors = FALSE)
    sortedOut <- merge(indexSpecDf, sortedOut, by = "fullName.y")
    # final sorted supertaxa list
    sortedOut$taxonID <- 0
    sortedOut$category <- "cat"
    sortedOut <- sortedOut[, c(
        "abbrName", "taxonID", "fullName.x", "species", "ncbiID",
        "sortedSupertaxon", "rank", "category")]
    colnames(sortedOut) <- c(
        "abbrName", "taxonID", "fullName", "ncbiID", "supertaxonID",
        "supertaxon", "rank", "category")
    sortedOut$ncbiID <- as.factor(sortedOut$ncbiID)
    sortedOut$supertaxon <- as.factor(sortedOut$supertaxon)
    sortedOut$category <- as.factor(sortedOut$category)
    return(sortedOut)
}


#' Reduce the filtered profile data into supertaxon level
#' @description Reduce data of the processed phylogenetic profiles from input
#' taxonomy rank into supertaxon level (e.g. from species to phylum)
#' @param filteredProfile dataframe contains the filtered profiles (see
#' ?parseInfoProfile, ?filterProfileData and ?filteredProfile)
#' @return A reduced dataframe contains only profile data for the selected
#' supertaxon rank. This dataframe contains only supertaxa and their value
#' (mVar1 & mVar2) for each gene.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export

reduceProfileCr <- function(filteredProfile) {
    if (is.null(filteredProfile)) stop("Profile data cannot be NULL!")

    # check if working with the lowest taxonomy rank; 1 for NO; 0 for YES
    flag <- 1
    if (length(unique(levels(as.factor(filteredProfile$numberSpec)))) == 1) {
        if (unique(levels(as.factor(filteredProfile$numberSpec))) == 1) {
            superDfExt <- filteredProfile[, c(
                "geneID", "supertaxon", "supertaxonID",
                "var1", "presSpec", "category", "orthoID", "var2", "paralog"
            )]
            flag <- 0
        }
    }
    if (flag == 1) {
        # get representative orthoID that has m VAR1 for each supertaxon
        mOrthoID <- filteredProfile[, c(
            "geneID", "supertaxon", "var1", "mVar1", "orthoID", "presSpec"
        )]
        mOrthoID <- subset(mOrthoID, mOrthoID$var1 == mOrthoID$mVar1)
        colnames(mOrthoID) <- c(
            "geneID", "supertaxon", "var1", "mVar1", "orthoID", "presSpec"
        )
        mOrthoID <- mOrthoID[!is.na(mOrthoID$orthoID), ]
        mOrthoID <- mOrthoID[, c("geneID", "supertaxon", "orthoID", "presSpec")]
        mOrthoID <- mOrthoID[!duplicated(mOrthoID[, seq_len(2)]), ]
        # get data set for PhyloProfile plotting (contains only supertaxa info)
        superDf <- subset(filteredProfile, select = c(
            "geneID", "supertaxon", "supertaxonID",
            "mVar1", "category", "mVar2", "paralog"
        ))
        superDf <- superDf[!duplicated(superDf), ]
        superDfExt <- merge(
            superDf, mOrthoID, by = c("geneID", "supertaxon"), all.x = TRUE
        )
        superDfExt <- superDfExt[, c(
            "geneID", "supertaxon", "supertaxonID",
            "mVar1", "presSpec", "category", "orthoID", "mVar2", "paralog"
        )]
        # rename mVar to var
        names(superDfExt)[names(superDfExt) == "mVar1"] <- "var1"
        names(superDfExt)[names(superDfExt) == "mVar2"] <- "var2"
    }
    return(superDfExt)
}
