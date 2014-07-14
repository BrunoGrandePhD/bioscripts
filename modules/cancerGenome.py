import MySQLdb
import re

class cancerGenomeDB():
    """handles connection and queries of cancer genome database instances"""
    #genome_fasta = None
    #for hg19 pass in this file instead: /projects/rmorin/common/genomes/hg19/GRCh37-lite.fa
    def __init__(self, database_name = None, database_host = None, database_user = None, database_password = None, genome_fasta = None):
        if genome_fasta:
            self.genome_fasta = genome_fasta
        else:
            self.genome_fasta = '/projects/rmorin/common/genomes/human_all.fasta'
            #genome_fasta = '/projects/rmorin/common/genomes/human_all.fasta'
        genome_fasta = self.genome_fasta
        if database_name:
            self.database = database_name
        else:
            self.database = 'lymphoma_meta_hg19'
        if database_password:
            self.password = database_password
        else:
            self.password = 'viewer'
        if database_user:
            self.user = database_user
        else:
            self.user = 'viewer'
        if database_host:
            self.host = database_host
        else:
            self.host = 'jango.bcgsc.ca'
        db = MySQLdb.connect(host=self.host, user=self.user, passwd=self.password, db=self.database)
        self.db = db

    def getSomaticSNVs(self,library=None):
        '''get all the somatic SNVs in the database'''
        cursor = self.db.cursor()
        query = "select somatic_mutation.id from somatic_mutation, event where event.id = somatic_mutation.event_id"

        if library:
            query = query + " and library_id = %s" % library.id
        #print query
        cursor.execute(query)
        snv_ids = []
        for result in cursor.fetchall():
            snv_ids.append(result[0])
        snvs = []
        for snv_id in snv_ids:
            snv = SNV(self.db,snv_id,genome_wide=True)
            snvs.append(snv)
        return snvs
    def getLibraries(self, type = None, sample_name = None, sample_id = None, subtype = None, sample_type = None, patient_id = None, treatment = None):
        """get libraries for all/some samples of a given (or all) library types such as genome, RNA-seq"""
        query = 'select library.id from library, sample where sample.id = library.sample_id '
        if type:
            if isinstance(type, list):
                query = query + ' and library_type in ('
                for i in type:
                    query = query + "'%s'," % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library_type = '%s'" % type
        if subtype:
            query = query + " and subtype ='%s'" % subtype
        if sample_type:
            query = query + " and sample_type = '%s'" % sample_type
        if patient_id:
            query = query + ' and patient_id = %s' % patient_id
        elif sample_name:
            query = query + ' and sample.sample_id = "%s"' % sample_name
        elif sample_id:
            query = query + ' and sample.id = %s' % sample_id
        elif treatment:
            query = query + " and treatment = '%s'" % treatment
        cursor = self.db.cursor()
        print query
        cursor.execute(query)
        libraries = []
        for lib_id in cursor.fetchall():
            library = Library(self.db, library_id=lib_id)
            print "LIBRARY:"
            print library
            libraries.append(library)

        return libraries
    def getSpliceSiteSNVs(self,library_id=None,gene=None):
        return(self._getSpliceSiteSNVs(library_id=library_id,gene=gene))
    def _getSpliceSiteSNVs(self,library_id=None,gene=None):
        cursor = self.db.cursor()
        query = "select splice_site_snv.id from splice_site_snv, event, gene_event where event.id = splice_site_snv.event_id and gene_event.event_id = event.id"
        if gene:
            query = "select splice_site_snv.id from splice_site_snv, event, gene_event where gene_event.event_id = event.id and event.id = splice_site_snv.event_id and gene_event.gene_id = %s" % (gene.id)
        if library_id:
            query = query + " and library_id = %s" % library_id

        cursor.execute(query)
        objects = []
        for result in cursor.fetchall():
            sss_id = result[0]
            #print "splice site ID: %s" % sss_id
            splice_obj = SpliceSiteSNV(self.db,sss_id)
            objects.append(splice_obj)
        return(objects)
    def getSpectrum(self,library_id):
        '''Determine the frequency of each of the six base changes in a library and return this result (ignoring upstream and downstream base context)'''
        #WARNING: CURRENTLY IGNORES SPLICE SITE SNVS!
        cursor = self.db.cursor()
        spectrum = {}
        CT = 0
        query = "select count(*) from mutation where library_id = %s and (ref_base ='C' and nref_base = 'T' or ref_base = 'G' and nref_base ='A')" % library_id
        #print query
        cursor.execute(query)
        try:
            CT = cursor.fetchone()[0]
        except TypeError:
            print "none for CT"
            pass
        query = "select count(*) from mutation where library_id = %s and (ref_base ='C' and nref_base = 'A' or ref_base = 'G' and nref_base ='T')" % library_id
        CA = 0
        cursor.execute(query)
        try:
            CA = cursor.fetchone()[0]
        except TypeError:
            print "none for CT"
            pass
        query = "select count(*) from mutation where library_id = %s and (ref_base ='C' and nref_base = 'G' or ref_base = 'G' and nref_base ='C')" % library_id
        CG = 0
        cursor.execute(query)
        try:
            CG = cursor.fetchone()[0]
        except TypeError:
            pass
        query = "select count(*) from mutation where library_id = %s and (ref_base ='T' and nref_base = 'A' or ref_base = 'A' and nref_base ='T')" % library_id
        TA = 0
        cursor.execute(query)
        try:
            TA = cursor.fetchone()[0]
        except TypeError:
            pass
        TG = 0
        query = "select count(*) from mutation where library_id = %s and (ref_base ='T' and nref_base = 'G' or ref_base = 'A' and nref_base ='C')" % library_id
        cursor.execute(query)
        try:
            TG = cursor.fetchone()[0]
        except TypeError:
            pass

        query = "select count(*) from mutation where library_id = %s and (ref_base ='T' and nref_base = 'C' or ref_base = 'A' and nref_base ='G')" % library_id
        cursor.execute(query)
        TC = 0
        try:
            TC = cursor.fetchone()[0]
        except TypeError:
            pass
        spectrum["CT"] = int(CT)
        spectrum["CA"] = int(CA)
        spectrum["CG"] = int(CG)
        spectrum["TA"] = int(TA)
        spectrum["TC"] = int(TC)
        spectrum["TG"] = int(TG)
        return spectrum
    def getFusions(self,library_id=None,confirmed=None):
        return(self._getFusions(library_id=library_id,confirmed=confirmed))
    def _getFusions(self,library_id=None,confirmed=None):
        '''Create fusion objects for all fusions (or those in a library)'''
        cursor = self.db.cursor()
        #print cursor
        query = "select fusion_transcript.id from fusion_transcript, event where event.id = fusion_transcript.event_id"

        if library_id:
            query = query + " and library_id = %s" % library_id
        if confirmed:
            query = query + " and confirmed = 'yes'"
        #print query
        cursor.execute(query)
        fusions = []
        for result in cursor.fetchall():
            #print "using: %s" % result[0]

            fusion = FusionTranscript(self.db,fusion_transcript_id=result[0])

            if fusion == None:
                print "error with %s" % result
                continue
            #fusion_rev = FusionTranscript(self.db,fusion_transcript_id=result[0],reverse_order=True)
            fusions.append(fusion)
            #fusions.append(fusion_rev)
        return fusions
    def getRearrangements(self, sample_name = None, to_validate=None, event_type=None, library=None):
        """get all rearrangements, can restrict to a given sample"""
        query = 'select distinct event_id from genomic_break, event, library, sample where sample.id = library.sample_id and library.id = event.library_id and event.id = genomic_break.event_id'
        if sample_name:
            query = query + " and sample.sample_id = '%s'" % sample_name
        elif library:
            query = query + " and library.id = %s" % library.id
        if to_validate:
            query = query + " and to_validate = 'yes'"
        if event_type:
            query = query + " and event.type = '%s'" % event_type
        query = query + " order by library_id, event_id"
        cursor = self.db.cursor()
        #print query
        #exit()
        cursor.execute(query)
        rearrangements = []
        for event_id in cursor.fetchall():
            rearrangement = Rearrangement(self.db, event_id[0])
            rearrangements.append(rearrangement)
        return rearrangements

    def getRearrangedGenes(self, sample_name = None, max_dist = 1):
        """get all genes with genomic breaks in or near them, as defined by max_dist"""
        query = "select chromosome, position from genomic_break, event, library, sample where sample.id = library.sample_id and library.id = event.library_id and genomic_break.event_id = event.id and sample.sample_id = '%s'" % sample_name
        cursor = self.db.cursor()
        cursor.execute(query)
        genes = []
        for region in cursor.fetchall():
            chr = region[0]
            pos = region[1]
            if max_dist:
                start = pos - max_dist
                end = pos + max_dist
                these_genes = self.getGenesInRegion(chr, start, end)
                #print these_genes
                genes += these_genes
        return genes

    def getGenesInRegion(self, chromosome, start, end, coding_region = None):
        """get all genes that are within the specified region"""
        cursor = self.db.cursor()
        query = "select id from gene where chromosome = '%s' and ((%s > start_position and %s < end_position) or (%s < start_position and %s > end_position) or (%s > start_position and %s < end_position))" % (chromosome,
         end,
         end,
         start,
         end,
         start,
         start)
        cursor.execute(query)
        gene_objects = []
        for result in cursor.fetchall():
            if coding_region:
                gene_id = result[0]
                query = "select count(*) from codon_gene where chromosome = '%s' and position >= %s and position <= %s" % (chromosome, start, end)
                cursor.execute(query)
                num_codons = cursor.fetchone()[0]
                if num_codons:
                    pass
                else:
                    continue
            gene_obj = Gene(self.db, gene_id=result[0])
            gene_objects.append(gene_obj)

        return gene_objects

    def getGenesWithinRegion(self, chromosome, start, end):
        '''Gets all genes entirely encompassed within a region
        (i.e., excluding those at the start and end).
        Method written by Bruno Grande.
        '''
        cursor = self.db.cursor()
        region = {
            'region_chromosome': chromosome,
            'region_start': start,
            'region_end': end
        }
        query = 'SELECT id FROM gene WHERE chromosome = "{region_chromosome}" AND start_position >= {region_start} AND end_position <= {region_end}'.format(**region)
        print query
        cursor.execute(query)
        gene_instances = []
        for gene_id in cursor.fetchall():
            gene_instances.append(Gene(self.db, gene_id=gene_id))
        return gene_instances

    def getGenesAtPosition(self, chromosome, position, padding=1000):
        '''Gets all overlapping genes at a given position.
        The padding keyword argument adds a margin around genes
        as a way to roughly include regulatory elements.
        Method written by Bruno Grande.
        '''
        cursor = self.db.cursor()
        locus = {
            'chromosome': chromosome,
            'position': position,
            'padding': padding
        }
        query = 'SELECT id FROM gene WHERE chromosome = "{chromosome}" AND start_position - {padding} < {position} AND end_position + {padding} > {position}'.format(**locus)
        print query
        cursor.execute(query)
        gene_instances = []
        for gene_id in cursor.fetchall():
            gene_instances.append(Gene(self.db, gene_id=gene_id))
        return gene_instances

    def getFusionGenes(self, sample_name = None, genome_confirmed = None):
        """get all genes in a CNV of a given maximum size"""
        query = 'select distinct gene_id from gene_event, event, fusion_transcript, library, sample where sample.id = library.sample_id and library.id = event.library_id and event.id = fusion_transcript.event_id and gene_event.event_id = event.id'
        query = query + " and sample.sample_id = '%s'" % sample_name
        if genome_confirmed:
            query = query + " and confirmed = 'yes'"
        cursor = self.db.cursor()
        cursor.execute(query)
        print query
        genes = []
        for id in cursor.fetchall():
            gene_obj = Gene(self.db, gene_id=id[0])
            genes.append(gene_obj)

        return genes

    def getSpliceSiteMutatedGenes(self, library_name = None, library_object = None, library_id = None):
        """get all genes in a CNV of a given maximum size"""
        query = 'select distinct gene_id from gene_event, event, splice_site_snv, library where library.id = event.library_id and event.id = splice_site_snv.event_id and gene_event.event_id = event.id'
        if library_name:
            if isinstance(library_name, list):
                query = query + ' and library_name in ('
                for i in library_name:
                    query = query + "'%s'," % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library_name = '%s'" % library_name
        elif library_id:
            if isinstance(library_id, list):
                query = query + ' and library.id in ('
                for i in library_id:
                    query = query + '%s,' % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library.id = '%s'"
        elif library_object:
            library = []
            if isinstance(library_object, list):
                for lib in library_object:
                    library.append(lib.id)

                print library
            else:
                library.append(library_object.id)
            return self.getSpliceSiteMutatedGenes(library_id=library)
        cursor = self.db.cursor()
        cursor.execute(query)
        genes = []
        for id in cursor.fetchall():
            gene_obj = Gene(self.db, gene_id=id[0])
            genes.append(gene_obj)

        return genes
    def getCNVs(self,max_size=2000000,limits=None,library_id=None):
        '''get all CNVs of a given maximum size'''
        cursor = self.db.cursor()
        query = "select cnv.id from cnv, event, library where event.id = cnv.event_id and cnv.type = 'somatic' and library.id = event.library_id and cnv.size <= %s and segment_state !=2" % max_size

        if library_id:
            query = query + " and library_id = %s" % library_id
        if limits:
            for key in limits:
                query = query + ' and %s %s' % (key, limits[key])
        print query
        #exit()
        cursor.execute(query)
        cnvs = []
        for result in cursor.fetchall():
            cnv_id = result[0]
            cnv_ob = CNV(self.db,cnv_id=cnv_id)
            cnvs.append(cnv_ob)
        return cnvs
    def getCNVGenes(self, library_name = None, library_object = None, library_id = None, max_size = None, state = None):
        """get all genes in a CNV of a given maximum size. """
        #WARNING, THIS IS BROKEN RIGHT NOW
        query = 'select distinct gene_id from gene_event, event, CNV, library where library.id = event.library_id and event.id = CNV.event_id and gene_event.event_id = event.id'
        if library_name:
            if isinstance(library_name, list):
                query = query + ' and library_name in ('
                for i in library_name:
                    query = query + "'%s'," % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library_name = '%s'" % library_name
        elif library_id:
            if isinstance(library_id, list):
                query = query + ' and library.id in ('
                for i in library_id:
                    query = query + '%s,' % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library.id = '%s'"
        elif library_object:
            library = []
            if isinstance(library_object, list):
                for lib in library_object:
                    library.append(lib.id)

                print library
            else:
                library.append(library_object.id)
            return self.getCNVGenes(library_id=library, max_size=max_size, state=state)
        if max_size:
            query = query + ' and size <= %i' % max_size
        if state:
            query = query + ' and segment_state in (%s)' % state
        else:
            query = query + ' and segment_state != 2'
        cursor = self.db.cursor()
        print query
        cursor.execute(query)
        genes = []
        for id in cursor.fetchall():
            gene_obj = Gene(self.db, gene_id=id[0])
            genes.append(gene_obj)

        return genes

    def getIndels(self, library_name = None, library_id = None, type = None, limits = None):
        return self._getIndels(library_name, library_id, type, limits)

    def _getIndels(self, library_name = None, library_id = None, type = None, limits = None):
        """get the indels from the database either from a given library or for all libraries with flexible limits allowed"""
        query = "select indel.id from library, indel, event where event.id = indel.event_id and event.library_id = library.id"
        cursor = self.db.cursor()
        if library_name:
            if isinstance(library_name, list):
                query = query + ' and library_name in ('
                for i in library_name:
                    query = query + "'%s'," % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library_name = '%s'" % library_name
        elif library_id:
            if isinstance(library_id, list):
                query = query + ' and library.id in ('
                for i in library_id:
                    query = query + '%s,' % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library.id = '%s'" % library_id
        else:
            query = 'select indel.id from indel, event where event.id = indel.event_id '
        if limits:
            for key in limits:
                query = query + ' and %s %s' % (key, limits[key])

        cursor.execute(query)
        indels = []
        for res in cursor.fetchall():
            indel = Indel(self.db, res[0])
            indels.append(indel)
        return indels

    def getMutations(self, library_name = None, library_id = None, type = None, limits = None):
        return self._getMutations(library_name, library_id, type, limits)

    def _getMutations(self, library_name = None, library_id = None, type = None, limits = None):
        """get the mutations from the database either from a given library or for all libraries with flexible limits allowed"""
        #query = 'select mutation.id from gene, mutation, library where mutation.gene = gene.ensembl_id and mutation.library_id = library.id and validation_outcome not in ("false","unclear","germline")'
        query = 'select mutation.id from gene, mutation, library where mutation.gene = gene.ensembl_id and mutation.library_id = library.id'
        #print query
        #print library_name
        cursor = self.db.cursor()
        if library_name:
            if isinstance(library_name, list):
                query = query + ' and library_name in ('
                for i in library_name:
                    query = query + "'%s'," % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library_name = '%s'" % library_name
        elif library_id:
            if isinstance(library_id, list):
                query = query + ' and library.id in ('
                for i in library_id:
                    query = query + '%s,' % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library.id = '%s'" % library_id
        else:
            query = 'select mutation.id from gene, mutation, library where mutation.library_id = library.id and mutation.gene = gene.ensembl_id'
            #print query
        if limits:
            for key in limits:
                query = query + ' and %s %s' % (key, limits[key])
        #query = 'select mutation.id from gene, mutation, library where mutation.library_id = library.id and mutation.gene = gene.ensembl_id'  #TEMPORARY
        #print query
        #exit()
        cursor.execute(query)
        snvs = []
        for res in cursor.fetchall():
            snv = SNV(self.db, res[0])
            #print snv
            snvs.append(snv)

        return snvs
    def getAlleleDetails(self,mutation_id=None,splice_site_snv_id=None):
        return(self._getAlleleDetails(mutation_id=mutation_id,splice_site_snv_id=splice_site_snv_id))
    def _getAlleleDetails(self,mutation_id=None,splice_site_snv_id=None):
        """Retrieve allele support information from database for a coding SNV or splice site SNV for all samples that have this information populated"""
        cursor = self.db.cursor()
        if mutation_id:
            query = "select sample_id, ref_count, nref_count, variant_allele_fraction, experiment_type, cellularity_estimate from allele_count where mutation_id = %s" % mutation_id
        elif splice_site_snv_id:
            query = "select sample_id, ref_count, nref_count, variant_allele_fraction, experiment_type, cellularity_estimate from allele_count,  splice_site_snv where splice_site_snv.event_id = allele_count.event_id and splice_site_snv.id = %s" % splice_site_snv_id
        else:
            print "error: you must specify an id for either a splice_site_snv or mutation"

        cursor.execute(query)
        sample_allele_data = {}
        for row in cursor.fetchall():
            (sample_id,ref_count,nref_count,vaf,exptype,cellularity) = row
            sample_id =int(sample_id)
            ref_count = int(ref_count)
            nref_count = int(nref_count)
            vaf = float(vaf)
            cellularity = float(cellularity)
            sample_allele_data[sample_id]={"ref_count":ref_count,"nref_count":nref_count,"variant_allele_fraction":vaf,"experiment_type":exptype,"cellularity_estimate":cellularity}
        return sample_allele_data
    def getMutatedGenes(self, only_validated = None, library_name = None, library_object = None, library_id = None, type = None, limits = None):
        """get mutated genes either in entire database or in a single, or list of libraries, either all mutations or coding only"""
        query = 'select distinct gene.id from gene, mutation, library where mutation.gene = gene.ensembl_id and mutation.library_id = library.id'

        if library_name:
            if isinstance(library_name, list):
                query = query + ' and library_name in ('
                for i in library_name:
                    query = query + "'%s'," % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library_name = '%s'" % library_name
        elif library_id:
            if isinstance(library_id, list):
                query = query + ' and library.id in ('
                for i in library_id:
                    query = query + '%s,' % i

                query = query.rstrip(',')
                query = query + ')'
            else:
                query = query + " and library.id = '%s'"
        else:
            if library_object:
                library = []
                if isinstance(library_object, list):
                    for lib in library_object:
                        library.append(lib.id)

                    print library
                else:
                    library.append(library_object.id)
                return self.getMutatedGenes(library_id=library, type=type, limits=limits, only_validated=only_validated)
            query = 'select distinct gene.id from gene, mutation where mutation.gene = gene.ensembl_id'
        if limits:
            for key in limits:
                query = query + ' and %s %s' % (key, limits[key])
        if only_validated:
            #filter on genes with validated_somatic = 'yes' in gene table
            query = query + "and validated_somatic = 'yes'"
        print query
        cursor = self.db.cursor()
        cursor.execute(query)
        genes = []
        for id in cursor.fetchall():
            gene_obj = Gene(self.db, gene_id=id[0])
            genes.append(gene_obj)

        return genes

    def addProtein(self, ensembl_protein_id, length, ensembl_transcript_id = None, transcript_id = None):
        """Add a new protein to the protein table (first check if it exists)"""
        cursor = self.db.cursor()
        check_query = "select id from protein where ensembl_id = '%s'" % ensembl_protein_id
        cursor.execute(check_query)
        results = cursor.fetchone()
        if results:
            print '%s already populated' % ensembl_protein_id
            return results[0]
        else:
            if ensembl_transcript_id:
                transcript = Transcript(self.db, ensembl_transcript_id=ensembl_transcript_id)
                transcript_id = transcript.id
            query = "insert into protein (transcript_id,ensembl_id,length) values(%s,'%s',%s)" % (transcript_id, ensembl_protein_id, length)
            cursor.execute(query)
            check_query = "select id from protein where ensembl_id = '%s'" % ensembl_protein_id
            cursor.execute(check_query)
            results = cursor.fetchone()
            return results[0]

    def addExon(self, ensembl_exon_id, transcript_start, transcript_end, genome_start, genome_end, strand, phase, end_phase, ensembl_transcript_id = None, transcript_id = None):
        """Add a new exon to the Exon table (first check if it exists).  Note, check both ensembl_id and transcript_id as some transcripts share the same exon!"""
        cursor = self.db.cursor()
        if ensembl_transcript_id:
            transcript = Transcript(self.db, ensembl_transcript_id=ensembl_transcript_id)
            print transcript
            transcript_id = transcript.id
        check_query = "select id from exon where ensembl_id = '%s' and transcript_id = %s" % (ensembl_exon_id, transcript_id)
        cursor.execute(check_query)
        results = cursor.fetchone()
        if results:
            print '%s already populated' % ensembl_exon_id
            return results[0]
        else:
            query = "insert into exon (ensembl_id,transcript_id,transcript_start,transcript_end,genome_start,genome_end,strand,phase,end_phase) values('%s',%s,%s,%s,%s,%s,'%s','%s','%s')" % (ensembl_exon_id,
             transcript_id,
             transcript_start,
             transcript_end,
             genome_start,
             genome_end,
             strand,
             phase,
             end_phase)
            cursor.execute(query)
            cursor.execute(check_query)
            results = cursor.fetchone()
            return results[0]

    def addProteinRegion(self, accession, source, db_id, description, protein_start, protein_end, protein_id = None, ensembl_protein_id = None):
        """Add annotations for a given protein from various sources"""
        cursor = self.db.cursor()
        if ensembl_protein_id:
            protein = Protein(self.db, ensembl_protein_id=ensembl_protein_id)
            protein_id = protein.id
        query = "insert into protein_region (protein_id,accession,source,db_id,description,protein_start,protein_end) values(%s,'%s','%s','%s','%s',%s,%s)" % (protein_id,
         accession,
         source,
         db_id,
         description,
         protein_start,
         protein_end)
        cursor.execute(query)
        cursor.execute('select last_insert_id()')
        region_id = cursor.fetchone()[0]
        return region_id

    def addTranscript(self, ensembl_gene_id, ensembl_transcript_id, length, cds_start, cds_end, refseq_id = None):
        """Add a new transcript to the Transcript table (first check if it exists)"""
        cursor = self.db.cursor()
        check_query = "select id from transcript where ensembl_id = '%s'" % ensembl_transcript_id
        cursor.execute(check_query)
        results = cursor.fetchone()
        if results:
            print '%s already populated' % ensembl_transcript_id
            return results[0]
        else:
            gene = Gene(self.db, ensembl_id=ensembl_gene_id)
            gene_id = gene.id
            query = "insert into transcript(gene_id,ensembl_id,length,cds_start,cds_end,refseq) values(%s,'%s',%i,%i,%i,'%s')" % (gene_id,
             ensembl_transcript_id,
             length,
             cds_start,
             cds_end,
             refseq_id)
            print query
            cursor.execute(query)
            cursor.execute(check_query)
            results = cursor.fetchone()
            has_id = results[0]
            return has_id

    def parse_snv(anno_line):
        """A function to parse the data out of the standard SNV format and return it as a dictionary"""
        data = {}
        anno_line = anno_line.rstrip('\n')
        anno_cols = anno_line.split(' ')
        chrpos = anno_cols[0]
        chr, pos = chrpos.split(':')
        data['chromosome'] = chr
        data['position'] = pos
        data['gene'] = anno_cols[1]
        ref_base = anno_cols[7]
        mut_base = anno_cols[8]
        anno_type = anno_cols[13]
        data['base_change'] = '%s>%s'
        data['base_change'] = data['base_change'] % (ref_base, mut_base)
        data['annotation'] = anno_cols[14]
        if anno_type == 'CODING':
            data['protein_altering'] = 'yes'
        else:
            data['protein_altering'] = 'no'
        return (data, anno_type)
    def parse_cnv(self,cnv_line,file_format):
        cnv_data = {}
        if file_format == "corbett":
            #23 1179669 1640633 3
            cnv_line = cnv_line.rstrip('\n')
            cnv_cols = cnv_line.split("\t")
            chr = cnv_cols[0]
            if chr == '23':
                chr = 'X'
            if chr == '24':
                chr = 'Y'
            chr = "chr%s" % chr
            print "chromosome: %s" % chr
            cnv_data['chromosome'] = chr
            cnv_data['start'] = int(cnv_cols[1])
            cnv_data['end'] = int(cnv_cols[2])
            if cnv_cols[3] == "HOMD":
                cnv_cols[3] = 0
            elif cnv_cols[3] == "HETD":
                cnv_cols[3] = 1
            elif cnv_cols[3] == "NEUT":
                cnv_cols[3] = 2
            elif cnv_cols[3] == "GAIN":
                cnv_cols[3] = 3
            elif cnv_cols[3] == "AMP":
                cnv_cols[3] = 4
            elif cnv_cols[3] == "HLAMP":
                cnv_cols[3] = 5
            cnv_data['state'] = int(cnv_cols[3])
        elif file_format == "dnacopy":
            #ID chrom loc.start loc.end num.mark seg.mean
            #some hard thresholds
            min_markers = 6
            gain_thresh = 0.5
            loss_thresh = -0.5
            cnv_line = cnv_line.rstrip()
            cnv_cols = cnv_line.split(" ")
            if cnv_cols[0] == "ID":
                return cnv_data

            print cnv_line
            chr = cnv_cols[2]
            chr = "chr" + chr
            if cnv_cols[5] == 0 or cnv_cols[6] == "NA":
                #not sure what causes this
                return cnv_data
            cnv_data['num_markers'] = int(cnv_cols[5])
            if cnv_data['num_markers']< min_markers:
                return cnv_data
            cnv_data['chromosome'] = chr
            cnv_data['start'] = int(float(cnv_cols[3]))
            cnv_data['end'] = int(float(cnv_cols[4]))
            cnv_data['segment_mean'] = float(cnv_cols[6])
            if cnv_data['segment_mean']  >= loss_thresh and cnv_data['segment_mean'] <= gain_thresh:
                cnv_data['state'] = 2
            elif cnv_data['segment_mean'] < loss_thresh:
                cnv_data['state'] = 1
            elif cnv_data['segment_mean'] > gain_thresh:
                cnv_data['state'] = 3
        elif file_format == "hmmcopy":
            #ID	chrom	start	end	num.mark	seg.mean	state.num	state.name
            #TENOM028001	1	1	2468000	2468	-0.027229	3	NEUT
            #TENOM028001	1	2468001	2512000	44	0.297031	4	GAIN
            cnv_line = cnv_line.rstrip('\n')
            cnv_cols = cnv_line.split("\t")
            if cnv_cols[0] == "ID":
                return cnv_data
            chr = cnv_cols[1]
            if chr == '23':
                chr = 'X'
            if chr == '24':
                chr = 'Y'
            chr = "chr%s" % chr
            print "chromosome: %s" % chr
            cnv_data['chromosome'] = chr
            cnv_data['start'] = int(float(cnv_cols[2]))
            cnv_data['end'] = int(float(cnv_cols[3]))
            cnv_data['state'] = int(cnv_cols[6])-1
        elif file_format == "titan":
            # Sample  Chromosome  Start_Position(bp)  End_Position(bp)    Length(bp)  Median_Ratio    \
            #       Median_logR TITAN_state TITAN_call  Copy_Number MinorCN MajorCN Clonal_Cluster  Clonal_Frequency
            # HY1972  1   13868   739528  725661  0.72    0.01    3   NLOH    2   0   2   4   0.58
            # HY1972  1   741267  829637  88371   0.67    0.14    9   ALOH    3   0   3   4   0.58
            cnv_cols = cnv_line.rstrip().split("\t")
            # Skip header line
            if cnv_cols[0] == "Sample":
                return cnv_data
            cnv_data['chromosome'] = cnv_cols[1]
            cnv_data['start'] = int(float(cnv_cols[2]))
            cnv_data['end'] = int(float(cnv_cols[3]))
            cnv_data['state'] = int(cnv_cols[9])
        else:
            print "error, no file type specified"
            exit()
        return(cnv_data)
    def addSomaticMutation(self, library_id, chromosome, position, base_change):
        """Load mutation from genome-wide mutation calls into database"""
        cursor = self.db.cursor()
        query = "select count(*) from somatic_mutation, event where event.id = somatic_mutation.event_id and library_id = %s and chromosome = '%s' and position = %s" % (library_id, chromosome, position)
        cursor.execute(query)
        num = cursor.fetchone()[0]
        if num:
            print 'already populated, doing nothing'
            return 0
        #add new event, then add new mutation
        query = "insert into event (library_id,type) values (%s,'SNV')" % library_id
        print query
        cursor.execute(query)
        query = "select last_insert_id()"
        cursor.execute(query)
        event_id = cursor.fetchone()[0]
        query = "insert into somatic_mutation (event_id,chromosome,position,base_change) values(%s,'%s','%s','%s')" % (event_id,
         chromosome,
         position,
         base_change)
        cursor.execute(query)
        print "added %s to event and %s:%s to somatic_mutation" % (event_id,chromosome,position)
    def addAlleleDetails(self, sample_id,chromosome,position,parent_library,ref_count,nref_count,mutation_type="coding",experiment_type=None,cellularity_estimate=None):
        """Add details for the two alleles for a given mutation to the database and link to Mutation and Sample tables.
        Note, data can refer to mutations not in the same library as the sample (e.g. in primary/met pairs). The library the mutation was called in is here referred to as the parent_library.
        The allele fraction data need not come from the bam file for library. For example, it could come from another sample from the same patient. This code will populate the allele_count table."""
        cursor = self.db.cursor()
        library_id = parent_library.id
        if mutation_type == "coding":
            mutation_query = "select id from mutation where library_id = %s and chromosome = '%s' and position = %s" % (library_id,chromosome,position)
            cursor.execute(mutation_query)
            mutation_id = cursor.fetchone()[0]
            check_query = "select count(*) from allele_count where allele_count.mutation_id = %s and sample_id = %s" % (mutation_id,sample_id)
        elif mutation_type == "splice":
            mutation_query = "select event_id from splice_site_snv, event where event.id = splice_site_snv.event_id and library_id = %s and chromosome = '%s' and position = %s" % (library_id,chromosome,position)
            cursor.execute(mutation_query)
            event_id = cursor.fetchone()[0]
            check_query = "select count(*) from allele_count where allele_count.event_id = %s and sample_id = %s" % (event_id,sample_id)
        #check that the details for this mutation have not already been added

        cursor.execute(check_query)
        num_records = cursor.fetchone()[0]
        if num_records >0:
            #this data is already in the database
            return()
        #insert data into allele_count table
        if not cellularity_estimate:
            cellularity_estimate = 0
        if not experiment_type:
            experiment_type = 'exome'
        variant_allele_fraction = float(nref_count)/(float(nref_count) + float(ref_count))
        #this code will be much cleaner when the mutation table gets linked through the event table
        if mutation_type == 'coding':
            insert_query = "insert into allele_count(mutation_id,sample_id,ref_count,nref_count,variant_allele_fraction,experiment_type,cellularity_estimate) values(%s,%s,%s,%s,%f,'%s',%f)" % (mutation_id,sample_id,ref_count,nref_count,variant_allele_fraction,experiment_type,cellularity_estimate)
        elif mutation_type == 'splice':
            insert_query = "insert into allele_count(event_id,sample_id,ref_count,nref_count,variant_allele_fraction,experiment_type,cellularity_estimate) values(%s,%s,%s,%s,%f,'%s',%f)" % (event_id,sample_id,ref_count,nref_count,variant_allele_fraction,experiment_type,cellularity_estimate)
        cursor.execute(insert_query)

    def addMutation(self, library_id, chromosome, position, ensembl_gene_id, base_change, status=None, protein_altering=None, annotation=None, validation_outcome=None, cdna_change = None, to_validate=None,identifiers=None, splice_site=None, mutation_seq_probability=None, triplet=None, tumour_ref = None, tumour_nref = None, normal_ref = None, normal_nref = None,transcript=None,sift_score=None,polyphen_score=None,mutation_ass_score=None,ref_base=None,nref_base=None):
        """Load new mutation into the database for a given library"""
        cursor = self.db.cursor()

        if splice_site:
            #work with different set of tables in the db
            gene_obj = Gene(self.db,ensembl_id=ensembl_gene_id)

            query = "select count(*) from splice_site_snv, event where splice_site_snv.event_id = event.id and event.library_id = %s and chromosome = '%s' and position = %s" % (library_id,chromosome,position)
            print query
            cursor.execute(query)
            num = cursor.fetchone()[0]
            if num:
                print 'already populated, doing nothing'
                return 0
            query = "insert into event (library_id,type) values('%s','SNV')" % library_id
            print query
            cursor.execute(query)
            query = "select last_insert_id()"
            cursor.execute(query)
            event_id = cursor.fetchone()[0]
            if mutation_seq_probability:
                query = "insert into splice_site_snv (chromosome, position, event_id, base_change,mutation_seq_prob,triplet_context,reference_base_count,nonreference_base_count,germline_reference_base_count,germline_nonreference_base_count) values('%s','%s','%s','%s',%f,'%s',%s,%s,%s,%s)" % (chromosome,
                position,
                event_id,
                base_change,
                mutation_seq_probability,
                triplet,
                tumour_ref,
                tumour_nref,
                normal_ref,
                normal_nref)
            else:
                query = "insert into splice_site_snv (chromosome, position, event_id, base_change,triplet_context,reference_base_count,nonreference_base_count,germline_reference_base_count,germline_nonreference_base_count) values('%s','%s','%s','%s','%s',%s,%s,%s,%s)" % (chromosome,
                position,
                event_id,
                base_change,
                triplet,
                tumour_ref,
                tumour_nref,
                normal_ref,
                normal_nref)
                cursor.execute(query)
            print query
            try:
                test = gene_obj.id
            except AttributeError:
                print "Gene not in db, inserting SNV with no gene"
                return 0
            query = "insert into gene_event(event_id,gene_id,nature,effect) values(%s,%s,'somatic','splice')" % (event_id,gene_obj.id)
            cursor.execute(query)
            print query
            print "inserted new event record (%s), splice_site_snv record and gene_event record pointing to %s" % (event_id,gene_obj.id)

        else:
            query = "select count(*) from mutation where library_id = %s and chromosome = '%s' and position = %s" % (library_id, chromosome, position)
            cursor.execute(query)
            num = cursor.fetchone()[0]
            if num:
                print 'already populated, doing nothing'
                return 0
            if validation_outcome == 'somatic' or validation_outcome == 'known_somatic':
                status = 'valid_somatic'
            if mutation_seq_probability:
                query = "insert into mutation (library_id,type,chromosome,position,gene,status,protein_altering,base_change,annotation,validation_outcome,to_validate,known_identifiers,cdna_change,mutation_seq_prob,triplet_context,reference_base_count,nonreference_base_count,germline_reference_base_count,germline_nonreference_base_count,transcripts,sift_score,polyphen_score,mutation_ass_score) values(%s,'snv','%s',%s,'%s','%s','%s','%s','%s','%s','%s','%s','%s',%f,'%s',%s,%s,%s,%s,'%s',%f,%f,%f)" % (library_id,
                chromosome,
                position,
                ensembl_gene_id,
                status,
                protein_altering,
                base_change,
                annotation,
                validation_outcome,
                to_validate,
                identifiers,
                cdna_change,
                mutation_seq_probability,
                triplet,
                tumour_ref,
                tumour_nref,
                normal_ref,
                normal_nref,
                transcript,
                sift_score,
                polyphen_score,
                mutation_ass_score)
            else:
                print chromosome
                print position
                print ensembl_gene_id
                print status
                print protein_altering
                print base_change
                print annotation
                print validation_outcome
                print to_validate
                print identifiers
                print cdna_change
                print triplet
                print tumour_ref
                print tumour_nref
                print normal_ref
                print normal_nref
                print transcript
                print sift_score
                print polyphen_score
                query = """insert into mutation (
                library_id,
                type,
                chromosome,
                position,
                gene,
                status,
                protein_altering,
                base_change,
                annotation,
                validation_outcome,
                to_validate,
                known_identifiers,
                cdna_change,
                triplet_context,
                reference_base_count,
                nonreference_base_count,
                germline_reference_base_count,
                germline_nonreference_base_count,
                transcripts,
                sift_score,
                polyphen_score)
                values(%s,'snv','%s',%s,'%s','%s','%s','%s','%s','%s','%s','%s','%s','%s',%s,%s,%s,%s,'%s',%f,%f)""" % (
                library_id,
                chromosome,
                position,
                ensembl_gene_id,
                status,
                protein_altering,
                base_change,
                annotation,
                validation_outcome,
                to_validate,
                identifiers,
                cdna_change,
                triplet,
                tumour_ref,
                tumour_nref,
                normal_ref,
                normal_nref,
                transcript,
                sift_score,
                polyphen_score)
            print "adding mutation_______________"
            print query
            cursor.execute(query)
            if nref_base:
                query = "select last_insert_id()"
                cursor.execute(query)
                mutation_id = cursor.fetchone()[0]
                query = "update mutation set nref_base = '%s', ref_base= '%s' where id = %s" % (nref_base,ref_base,mutation_id)
                cursor.execute(query)

            return 1
    def addLohSegment(self, library_id, chromosome, segment_start, segment_end, size, bin_count, copy_number, loh_state, Major_allele_count, minor_allele_count):
        """Load LOH segment from Apolloh into the database for a given library"""
        cursor = self.db.cursor()
        query = "select count(*) from loh, event where event.id = loh.event_id and library_id = %s and chromosome = '%s' and segment_start = %s" % (library_id, chromosome, segment_start)
        cursor.execute(query)
        num = cursor.fetchone()[0]
        if num:
            print 'already populated, doing nothing'
            return 0
        #first insert a new loh event, get the event_id, then load details into the table
        query = "insert into event (library_id,type) values(%s,'loh')" % library_id
        cursor.execute(query)
        query = "select last_insert_id()"
        cursor.execute(query)
        event_id = cursor.fetchone()[0]

        query = "insert into loh (event_id, chromosome, segment_start, segment_end, size, bin_count, copy_number, loh_state, Major_allele_count, minor_allele_count) values (%s,'%s',%s,%s,%s,%s,%s,'%s',%s,%s)" % (event_id,chromosome, segment_start, segment_end, size, bin_count, copy_number, loh_state, Major_allele_count, minor_allele_count)
        cursor.execute(query)
        return 1
    def loadTranscripts(self,filepath):
        '''Loads Ensembl transcripts into database while properly linking
        to the corresponding entry in the `gene` table.
        Method written by Bruno Grande.
        '''
        # Example row in the file
        # ENST00000593546 ENSG00000268612         2255    32      850
        cursor = self.db.cursor()
        file_handle = open(filepath)
        column_names = ('transcript_ens_id', 'gene_ens_id', 'length', 'cds_start', 'cds_end')
        for line in file_handle:
            transcript_dict = dict(zip(column_names, line.rstrip().split('\t')))
            # Obtain gene ID from the gene table using the Ensembl ID
            query = 'SELECT id FROM gene WHERE ensembl_id = "{gene_ens_id}"'.format(**transcript_dict)
            print query
            cursor.execute(query)
            transcript_dict['gene_id'] = int(cursor.fetchone()[0])
            # Add the transcript to the database while linking to the gene
            query = 'INSERT INTO transcript (ensembl_id, gene_id, length, cds_start, cds_end) VALUES ("{transcript_ens_id}", {gene_id}, {length}, {cds_start}, {cds_end})'.format(**transcript_dict)
            print query
            cursor.execute(query)

    def loadExons(self,file):
        '''load all ensembl exons into the db and link to transcript table'''
        #ENSE00003051371	ENST00000593546	1	272	26597180	26597451	1	-1	1
        cursor = self.db.cursor()
        handle = open(file,"r")
        for line in handle:
            line = line.rstrip()
            (ense,enst,tstart,tend,gstart,gend,strand,phase,end_phase) = line.split()
            query = "select id from transcript where ensembl_id = '%s'" % (enst)
            cursor.execute(query)
            trans_id = int(cursor.fetchone()[0])
            query = "insert into exon (transcript_id,ensembl_id,transcript_start,transcript_end,genome_start,genome_end,strand,phase,end_phase) values(%s,'%s',%s,%s,%s,%s,'%s','%s','%s')" % (trans_id,ense,tstart,tend,gstart,gend,strand,phase,end_phase)
            print query
            cursor.execute(query)
    def loadProteins(self,file):
        '''Load all ensembl proteins into the db and link to transcript table'''
        #ENSP00000468826 ENST00000593546 272
        #ensembl_id transcript length
        cursor = self.db.cursor()
        handle = open(file,"r")
        for line in handle:
            line = line.rstrip()
            (ensp,enst,length) = line.split()
            query = "select id from transcript where ensembl_id = '%s'" % (enst)
            print query
            cursor.execute(query)
            trans_id = int(cursor.fetchone()[0])
            query = "insert into protein (transcript_id,ensembl_id,length) values(%s,'%s',%s)" % (trans_id,ensp,length)
            print query
            cursor.execute(query)
    def loadCNVs(self,cnv_file,sample_name=None,library_name=None,file_format="hmmcopy"):
        if sample_name:
            lib_obj = self.getLibraries(type='genome',sample_name=sample_name,sample_type='tumour')[0]
        elif library_name:
            lib_obj = Library(self.db,library_name=library_name)
        else:
            print "error, no library or sample information given"
            exit()
        print lib_obj
        cursor = self.db.cursor()
        #check whether this sample has CNVs already
        query = "select sample_id from library where library.id = %i" % lib_obj.id
        cursor.execute(query)
        sid = cursor.fetchone()[0]
        query = "select count(*) from event, library where library.id = event.library_id and sample_id =  %i and type = 'CNV'" % sid
        cursor.execute(query)
        num = cursor.fetchone()
        if num[0] > 0:
            print "CNVs already populated for this library, skipping"
            return
        cnv_data = open(cnv_file,"r")
        cnv_lines = cnv_data.readlines()
        library_id = lib_obj.id
        counter = 0
        for line in cnv_lines:
            cnv_details = self.parse_cnv(line,file_format)
            if not cnv_details.has_key("chromosome"):
                continue
            print cnv_details
            query = "select count(*) from event, cnv where cnv.event_id = event.id and chromosome = '%s' and segment_start = %i and library_id = %i" % (cnv_details['chromosome'],cnv_details['start'],library_id)
            cursor.execute(query)
            result= cursor.fetchone()[0]
            if result > 0:
                print "already in database"
                return
            #first add the "event" to the event table (type=cnv,library_id=library_id)
            query = "insert into event (library_id,type) VALUES('%s','%s')" % (library_id,"cnv")
            print query
            print cnv_details
            cursor.execute(query)
            query = "select max(id) from event where library_id = %i" % library_id
            cursor.execute(query)
            event_id = cursor.fetchone()[0]
            print "event: %i" % event_id
            #now add cnv details
            size = cnv_details['end'] - cnv_details['start'] +1
            query = "insert into cnv (event_id,chromosome,segment_start,segment_end,segment_state,size,type) VALUES(%i,'%s',%i,%i,%i,%i,'%s')" % (event_id,cnv_details['chromosome'],cnv_details['start'],cnv_details['end'],cnv_details['state'],size,'somatic')
            cursor.execute(query)
            counter+=1
            #now populate 'gene_event' table based on genes in region
            genes = self.getGenesInRegion(cnv_details['chromosome'],cnv_details['start'],cnv_details['end'])
            effect = ""
            if cnv_details['state'] > 2:
                effect = "gain"
            elif cnv_details['state'] < 2:
                effect = "loss"
            if effect != "":
                if len(genes) == 0:
                    continue
                for gene in genes:
                    query = "insert into gene_event (event_id,gene_id,effect) values(%i,%i,'%s')"
                    query = query % (event_id,gene.id,effect)
                    print query
                    cursor.execute(query)
        print "loaded %s events" % counter
    def loadIndels(self, indel_file, file_format, tumour_library, analysis_type):
        """Read transabyss output files and import into database, analysis_type is either bubble_pop or split_alignment, file_format is either validator or full"""
        #Warning. This is now broken and needs to be replaced to work with new addIndel function
        try:
            indel_handle = open(indel_file, 'r')
        except IOError:
            print "error, can't open file"
            exit()

        if file_format == 'validator':
            for line in indel_handle:
                line = line.rstrip('\n')
                indel_details = {}
                vals = line.split()
                if vals[0] == 'rearrangement':
                    continue
                indel_details['event_type'] = vals[0]
                m = re.search('(chr.+):([0-9]+)-([0-9]+)', vals[1])
                chrom = m.group(1)
                start = m.group(2)
                end = m.group(3)
                indel_details['chromosome'] = chrom
                indel_details['start'] = start
                indel_details['end'] = end
                indel_details['length'] = vals[2]
                indel_details['contig'] = vals[3]
                indel_details['read_support_tumour'] = vals[4]
                indel_details['read_support_normal'] = vals[5]
                indel_details['validator_result'] = vals[7]
                if vals[7] == 'somatic':
                    self.addIndel(indel_details, tumour_library.id, analysis_type, validator=True)

    def addCoverage(self,library,coverage_type,total_bases,region_size=None,gene=None,cutoff=None,fraction_above_cutoff=None,chromosome=None,position=None):
        '''Add details of the sequence coverage for a gene, coordinate or library'''
        cursor = self.db.cursor()

        if coverage_type == "library":
            #will record the average coverage across all exonic regions
            #for speed, this can just be parsed from the end of the summarized coverage file (output by summarize_bedcoverage.pl)
            average_depth = total_bases / region_size
            library_id = library.id
            if cutoff:
                query = "insert into coverage (total_bases,average_depth,cutoff,fraction_above_cutoff) values(%f,%s,%f,%f)" % (total_bases,average_depth,cutoff,fraction_above_cutoff)
            else:
                query = "insert into coverage (total_bases,average_depth) values(%f,%s)" % (total_bases,average_depth)
            cursor.execute(query)
            query ="select last_insert_id()"
            cursor.execute(query)
            cov_id = cursor.fetchone()[0]
            query = "insert into library_coverage (library_id,coverage_id) values(%s,%s)" % (library_id,cov_id)
            cursor.execute(query)
        elif coverage_type == "gene":
            average_depth = total_bases / region_size
            gene_id = gene.id
            library_id = library.id
            #check if already present
            query = "select count(*) from coverage, gene_coverage where coverage.id = gene_coverage.coverage_id and gene_id = %s and library_id = %s" % (gene_id,library_id)
            cursor.execute(query)

            num = cursor.fetchone()[0]
            if num > 0:
                print "error, this is already in the db, skipping"
                return()

            #will record the average coverage across all exons of a specified gene
            #first add coverage information
            if cutoff:
                query = "insert into coverage (total_bases,average_depth,cutoff,fraction_above_cutoff) values(%f,%s,%f,%f)" % (total_bases,average_depth,cutoff,fraction_above_cutoff)
            else:
                query = "insert into coverage (total_bases,average_depth) values(%f,%s)" % (total_bases,average_depth)
            cursor.execute(query)
            query ="select last_insert_id()"
            cursor.execute(query)
            cov_id = cursor.fetchone()[0]
            query = "insert into gene_coverage (gene_id,library_id,coverage_id) values(%s,%s,%s)" % (gene_id,library_id,cov_id)
            cursor.execute(query)
        elif coverage_type == "coordinate":
            library_id = library.id
            #will record the coverage at specified genomic coordinates only
            #check if already present
            query = "select count(*) from  position_coverage, position where position.id = position_coverage.position_id and chromosome = '%s' and position = %s and library_id = %s" % (chromosome, position,library_id)
            cursor.execute(query)

            num = cursor.fetchone()[0]
            if num > 0:
                print "error, this is already in the db, skipping"
                return()

            average_depth = total_bases / 1

            #first get position id from position table or create one
            query = "select id from position where chromosome = '%s' and position = %s" % (chromosome, position)
            cursor.execute(query)
            try:
                pos_id = cursor.fetchone()[0]
            except TypeError:
                #need to add this position
                query = "insert into position (chromosome, position) values ('%s',%s)" % (chromosome,position)
                cursor.execute(query)
                query ="select last_insert_id()"
                cursor.execute(query)
                pos_id = cursor.fetchone()[0]
            if cutoff:
                query = "insert into coverage (total_bases,average_depth,cutoff,fraction_above_cutoff) values(%f,%s)" % (total_bases,average_depth,cutoff,fraction_above_cutoff)
            else:
                query = "insert into coverage (total_bases,average_depth) values(%f,%s)" % (total_bases,average_depth)
            cursor.execute(query)
            query ="select last_insert_id()"
            cursor.execute(query)
            cov_id = cursor.fetchone()[0]

            query = "insert into position_coverage (position_id, library_id, coverage_id) values (%s,%s,%s)" % (pos_id,library_id,cov_id)
            print query
            cursor.execute(query)

    def addIndel(self,library_id,chromosome,start,end,ref,alt,effect,annotation,ensembl_id=None):
        cursor = self.db.cursor()
        query = "select indel.id from indel, event where event.id = indel.event_id and chromosome = '%s' and (start = %s or end = %s) and library_id = %s" % (chromosome,start,end,library_id)
        print query
        cursor.execute(query)
        try:
            db_id = cursor.fetchone()[0]
        except TypeError:
            pass
        else:
            print "skipping duplicate indel"
            return()
        #add event record then indel and finally gene_event details
        query = "insert into event (library_id,type) values('%s','indel')" % (library_id)
        cursor.execute(query)
        query = "select last_insert_id()"
        cursor.execute(query)
        event_id = cursor.fetchone()[0]
        query = "insert into indel (event_id,chromosome,start,end,ref,alt,annotation,effect) values(%s,'%s',%s,%s,'%s','%s','%s','%s')" % (event_id,chromosome,start,end,ref,alt,annotation,effect)
        cursor.execute(query)
        query = "select last_insert_id()"
        cursor.execute(query)
        indel_id = cursor.fetchone()[0]
        #add gene details
        if ensembl_id:
            try:
                gene_obj = Gene(self.db,ensembl_id=ensembl_id)
                id = gene_obj.id
            except AttributeError:
                print "error getting gene info for %s" % ensembl_id
                exit()
            query = "insert into gene_event (event_id,gene_id,nature) values (%s,%s,'somatic')" % (event_id,id)
            cursor.execute(query)
        return indel_id

    def addIndelLegacy(self, indel_details, library_id, indel_source, validator = None):
        cursor = self.db.cursor()
        query = "select indel.id, validator_result from indel, event where event.id = indel.event_id and chromosome = '%s' and (start = %s or end = %s) and library_id = %s" % (indel_details['chromosome'],
         indel_details['start'],
         indel_details['end'],
         library_id)
        cursor.execute(query)
        try:
            db_id, val_res = cursor.fetchone()
        except TypeError:
            pass
        else:
            print 'indel %s:%s-%s is a duplicate (%s,%s), skipping' % (indel_details['chromosome'],
             indel_details['start'],
             indel_details['end'],
             db_id,
             val_res)
            return ()

        if validator:
            genes = self.getGenesInRegion(indel_details['chromosome'], indel_details['start'], indel_details['end'], coding_region=True)
            if len(genes) == 0:
                return
            query = "insert into event (library_id,type) values(%i,'indel')" % library_id
            cursor.execute(query)
            query = 'select last_insert_id()'
            cursor.execute(query)
            event_id = cursor.fetchone()[0]
            query = "insert into indel (event_id,event_type,indel_source,chromosome,start,end,contig,length,read_support_tumour,read_support_normal,validator_result) values (%i,'%s','%s','%s',%s,%s,'%s',%s,%s,%s,'%s')" % (event_id,
             indel_details['event_type'],
             indel_source,
             indel_details['chromosome'],
             indel_details['start'],
             indel_details['end'],
             indel_details['contig'],
             indel_details['length'],
             indel_details['read_support_tumour'],
             indel_details['read_support_normal'],
             indel_details['validator_result'])
            print query
            cursor.execute(query)
            for gene in genes:
                query = "insert into gene_event (event_id,gene_id,detailed_annotation) values(%i,%i,'%s')"
                print query
                try:
                    query = query % (event_id, gene.id, indel_details['event_type'])
                except TypeError:
                    print 'eror with'
                    print gene.ensembl_id
                    continue

                print query
                cursor.execute(query)

        else:
            print 'this code is broken right now'
            exit()
            query = "insert into event (library_id,type) values(%i,'indel')" % library_id
            cursor.execute(query)
            query = 'select last_insert_id()'
            cursor.execute(query)
            event_id = cursor.fetchone()[0]
            try:
                query = "insert into indel (event_id, event_type, chromosome, start, end, contig, contig_len, contig_start, contig_end, length, ref, alt, event_reads, contig_reads, genome_reads, repeat_length, from_end, confirm_contig_region, within_simple_repeats, repeatmasker, within_segdup, one_read_opposite, indel_source) VALUES(%i,'%s','%s',%i,%i,'%s',%i,%i,%i,%i,'%s','%s',%i,%i,%i,%i,%i,'%s','%s','%s','%s','%s','%s')" % (event_id,
                 indel_details['event_type'],
                 indel_details['chromosome'],
                 indel_details['start'],
                 indel_details['end'],
                 indel_details['contig'],
                 indel_details['contig_len'],
                 indel_details['contig_start'],
                 indel_details['contig_end'],
                 indel_details['length'],
                 indel_details['ref'],
                 indel_details['alt'],
                 indel_details['event_reads'],
                 indel_details['contig_reads'],
                 indel_details['genome_reads'],
                 indel_details['repeat_length'],
                 indel_details['from_end'],
                 indel_details['confirm_contig_region'],
                 indel_details['within_simple_repeats'],
                 indel_details['repeatmasker'],
                 indel_details['within_segdup'],
                 indel_details['one_read_opposite'],
                 indel_source)
            except:
                print "one of the following values doesn't match the expected type:" + event_id + ' ' + indel_details['event_type'] + ' ' + indel_details['chromosome'] + ' ' + indel_details['start'] + ' ' + indel_details['end'] + ' ' + indel_details['contig'] + ' ' + indel_details['contig_len'] + ' ' + indel_details['contig_start'] + ' ' + indel_details['contig_end'] + ' ' + indel_details['length'] + ' ' + indel_details['ref'] + ' ' + indel_details['alt'] + ' ' + indel_details['event_reads'] + ' ' + indel_details['contig_reads'] + ' ' + indel_details['genome_reads'] + ' ' + indel_details['repeat_length'] + ' ' + indel_details['from_end'] + ' ' + indel_details['confirm_contig_region'] + ' ' + indel_details['within_simple_repeats'] + ' ' + indel_details['repeatmasker'] + ' ' + indel_details['within_segdup'] + ' ' + indel_details['one_read_opposite'] + ' ' + indel_source
                return 0

            cursor.execute(query)
            genes = get_genes_in_event(event_id, 'indel', indel_details['chromosome'], indel_details['start'], indel_details['end'])
            if len(genes) == 0:
                return
            for gene in genes:
                query = "insert into gene_event (event_id,gene_id,detailed_annotation) values(%i,%i,'%s')"
                try:
                    query = query % (event_id, gene[0], indel_details['event_type'])
                except TypeError:
                    print 'eror with'
                    print gene[0]
                    continue

                print query
                cursor.execute(query)

    def addGenomicBreak(self, library_id, chromosome, position, side, nature):
        '''Creates a new genomic_break entry in database for a given library.
        Returns the id for the newly created genomic_break entry.
        If the entry already exists, returns its ID.
        Method written by Bruno Grande.
        '''
        cursor = self.db.cursor()
        genomic_break = {
            'library_id': library_id,
            'chromosome': chromosome,
            'position': position,
            'side': side
        }

        # Checks if the genomic break already exists for this library
        query = 'SELECT genomic_break.id FROM genomic_break, event WHERE genomic_break.event_id = event.id AND event.library_id = "{library_id}" AND chromosome = "{chromosome}" AND position = {position} AND side = "{side}"'.format(**genomic_break)
        print query
        count = cursor.execute(query)
        if count == 1:
            print 'Genomic break already exists. Adding nothing.'
            genomic_break['id'] = cursor.fetchone()[0]
            return genomic_break['id']
        if count > 1:
            raise ValueError('Multiple genomic_break entries already exist with those exact '
                             'attributes. It\'s impossible to determine which one to select.')

        # If the genomic break doesn't already exist, create an event entry
        query = 'INSERT INTO event (library_id, type) VALUES ("{library_id}", "SV")'.format(**genomic_break)
        print query
        cursor.execute(query)
        query = 'SELECT LAST_INSERT_ID()'
        print query
        cursor.execute(query)
        genomic_break['event_id'] = cursor.fetchone()[0]

        # Now that we have created an event entry we can create a gene_event entry to find
        # possible affected genes
        genes = self.getGenesAtPosition(chromosome, position)
        for gene in genes:
            gene_id = gene.gene_id
            # CREATE GENE EVENT ENTRY!!
            query = 'INSERT INTO gene_event (event_id, gene_id, nature, effect) VALUES ({event_id}, {gene_id}, "{nature}", "disruption")'.format(event_id=genomic_break['event_id'], gene_id=gene_id, nature=nature)
            print query
            cursor.execute(query)

        # Now that we have the newly created event_id, we can create a genomic_break entry
        query = 'INSERT INTO genomic_break (chromosome, position, event_id, side) VALUES ("{chromosome}", {position}, {event_id}, "{side}")'.format(**genomic_break)
        print query
        cursor.execute(query)
        query = 'SELECT LAST_INSERT_ID()'
        print query
        cursor.execute(query)
        genomic_break['id'] = cursor.fetchone()[0]
        return genomic_break['id']

    def addSvCnv(self, library_id, chromosome, segment_start, segment_end, segment_state,
                 nature):
        '''Creates a new cnv entry in the database for a given library.
        This is for CNVs found by ABySS.
        Returns the id for the newly created cnv entry.
        If the entry already exists, returns its ID.
        Method written by Bruno Grande and Marija Jovanovic.
        '''
        cursor = self.db.cursor()
        cnv = {
            'library_id': library_id,
            'chromosome': chromosome,
            'segment_start': segment_start,
            'segment_end': segment_end,
            'segment_state': segment_state,
            'size': int(segment_end) - int(segment_start) + 1,
            'nature': nature
        }

        #Checks if the cnv already exists for this library
        query = 'SELECT cnv.id FROM cnv, event WHERE cnv.event_id = event.id AND event.library_id = "{library_id}" AND chromosome = "{chromosome}" AND segment_start = "{segment_start}" AND segment_end = "{segment_end}" AND segment_state = "{segment_state}" AND cnv.type = "{nature}"'.format(**cnv)
        print query
        count = cursor.execute(query)
        if count == 1:
            print 'BRUNO THIS ALREADY EXISTS. Adding nothing.'
            cnv['id'] = cursor.fetchone()[0]
            return cnv['id']
        if count > 1:
            raise ValueError('Multiple CNV enteries already exist with those exact'
                             'attributes. It\'s impossible to determine which one to select.')

        #If the CNV doesn't exist, create an event entry
        query = 'INSERT INTO event (library_id, type) VALUES ("{library_id}", "SV")'.format(**cnv)
        print query
        cursor.execute(query)
        query = 'SELECT LAST_INSERT_ID()'
        print query
        cursor.execute(query)
        cnv['event_id'] = cursor.fetchone()[0]

        # Now that we have created an event entry we can create a gene_event entry to find
        # possible affected genes
        genes = self.getGenesWithinRegion(chromosome, segment_start, segment_end)
        for gene in genes:
            gene_id = gene.gene_id
            #Find effect from segment state
            if cnv['segment_state'] == 3:
                effect = 'gain'
            if cnv['segment_state'] == 1:
                effect = 'loss'
            # CREATE GENE EVENT ENTRY!!
            query = 'INSERT INTO gene_event (event_id, gene_id, nature, effect) VALUES ({event_id}, {gene_id}, "{nature}", "{effect}")'.format(event_id=cnv['event_id'], gene_id=gene_id, nature=nature, effect=effect)
            print query
            cursor.execute(query)

        # Now that we have the newly created event_id, we can create a cnv entry.... YAY
        query = 'INSERT INTO cnv (event_id, chromosome, segment_start, segment_end, size, type, segment_state) VALUES ("{event_id}", "{chromosome}", "{segment_start}", "{segment_end}", "{size}", "{nature}", "{segment_state}")'.format(**cnv)
        print query
        cursor.execute(query)
        query = 'SELECT LAST_INSERT_ID()'
        print query
        cursor.execute(query)
        cnv['id'] = cursor.fetchone()[0]
        return cnv['id']

    def addStructuralVariant(self, library_id, break1_chromosome, break1_position, break1_side,
                             break2_chromosome, break2_position, break2_side, sv_type,
                             num_read_pairs, num_spanning_reads, status):
        '''Creates a new structural_variant entry in database combining two genomic_break entries.
        Returns the id for the newly created structural_variant entry.
        If the entry already exists, returns its ID.
        Method written by Bruno Grande.
        '''
        cursor = self.db.cursor()
        structural_variant = {
            'break1_chromosome': break1_chromosome,
            'break1_position': break1_position,
            'break1_side': break1_side,
            'break2_chromosome': break2_chromosome,
            'break2_position': break2_position,
            'break2_side': break2_side,
            'sv_type': sv_type,
            'num_read_pairs': num_read_pairs,
            'num_spanning_reads': num_spanning_reads,
            'status': status,
            'cnv_id': 'NULL'
        }

        # Figure out status (somatic, germline, or unknown)
        if structural_variant['status'] == 'other':
                nature = 'germline'
        elif structural_variant['status'] == 'notfound':
            nature = 'unknown'
        elif (structural_variant['status'] == 'somatic' or
              structural_variant['status'] == 'related'):
            nature = 'somatic'
        else:
            raise ValueError('Unrecognized SV status '
                             '(expected: other, notfound, somatic, or related)')
        structural_variant['nature'] = nature

        # Creates and/or obtains IDs for associated genomic_break entries
        structural_variant['break1_id'] = self.addGenomicBreak(library_id, break1_chromosome,
                                                               break1_position, break1_side,
                                                               nature)
        structural_variant['break2_id'] = self.addGenomicBreak(library_id, break2_chromosome,
                                                               break2_position, break2_side,
                                                               nature)

        # Creates and/or obtains IDs for associated cnv entries
        if (structural_variant['sv_type'] == 'duplication' or
                structural_variant['sv_type'] == 'deletion'):
            segment_state = 3 if structural_variant['sv_type'] == 'duplication' else 1
            structural_variant['cnv_id'] = self.addSvCnv(library_id, break1_chromosome,
                                                         break1_position, break2_position,
                                                         segment_state, nature)

        # Checks if the structural variant already exists
        query = 'SELECT id FROM structural_variant WHERE break1_id = {break1_id} AND break2_id = {break2_id} AND type = "{sv_type}"'.format(**structural_variant)
        print query
        count = cursor.execute(query)
        if count == 1:
            print 'Structural variant already exists. Nothing added.'
            structural_variant['id'] = cursor.fetchone()[0]
            return structural_variant['id']
        if count > 1:
            raise ValueError('Multiple structural_variant entries already exist with those exact '
                             'attributes. It\'s impossible to determine which one to select.')

        # If the structural variant doesn't already exist, create an structural_variant entry
        query = 'INSERT INTO structural_variant (type, break1_id, break2_id, cnv_id, num_read_pairs, num_spanning_reads, status) VALUES ("{sv_type}", {break1_id}, {break2_id}, {cnv_id}, {num_read_pairs}, {num_spanning_reads}, "{status}")'.format(**structural_variant)
        print query
        cursor.execute(query)
        query = 'SELECT LAST_INSERT_ID()'
        print query
        cursor.execute(query)
        structural_variant['id'] = cursor.fetchone()[0]
        return structural_variant['id']

    def addLibrary(self, library_name, library_type, sample_id):
        """Add new libraries to the database with minimal details"""
        cursor = self.db.cursor()
        check_query = "select count(*) from library where library_name = '%s'" % library_name
        cursor.execute(check_query)
        num = cursor.fetchone()[0]
        if num:
            query = "select id from library where library_name = '%s'" % library_name
            cursor.execute(query)
            id = cursor.fetchone()[0]
            return id
        query = "insert into library (library_name,library_type,sample_id) values('%s','%s',%s)" % (library_name, library_type, sample_id)
        print query
        cursor.execute(query)
        query = "select id from library where library_name = '%s'" % library_name
        cursor.execute(query)
        id = cursor.fetchone()[0]
        return id

    def addSample(self, sample_name, patient_id, sample_type):
        """Add a new sample for a patient"""
        cursor = self.db.cursor()
        check_query = "select count(*) from sample where sample_id = '%s'" % sample_name
        cursor.execute(check_query)
        num = cursor.fetchone()[0]
        if num:

            query = "select id from sample where sample_id = '%s' and sample_type = '%s'" % (sample_name, sample_type)
            print query
            cursor.execute(query)
            id = cursor.fetchone()[0]
            return id
        query = "insert into sample (sample_id,patient_id,sample_type) values ('%s',%s,'%s')" % (sample_name, patient_id, sample_type)
        print query
        cursor.execute(query)
        query = "select id from sample where sample_id = '%s' and sample_type = '%s'" % (sample_name, sample_type)
        print query
        cursor.execute(query)
        id = cursor.fetchone()[0]
        return id

    def addPatient(self, external_id):
        """Add new patient to the database, return id"""
        cursor = self.db.cursor()
        check_query = "select id from patient where res_id = '%s'" % external_id
        cursor.execute(check_query)
        res = cursor.fetchone()
        if res:
            pass
        else:
            query = "insert into patient (res_id) values ('%s')" % external_id
            cursor.execute(query)
        cursor.execute(check_query)
        id = cursor.fetchone()[0]
        return id

    def get_gene_types(self):
        gene_types = {}
        gene_list = []
        gene_types['EZH2'] = 'Histone_modifier'
        gene_list.append('EZH2')
        gene_types['CREBBP'] = 'Histone_modifier'
        gene_list.append('CREBBP')
        gene_types['MEF2B'] = 'Histone_modifier'
        gene_list.append('MEF2B')
        gene_types['MLL2'] = 'Histone_modifier'
        gene_list.append('MLL2')
        gene_types['EP300'] = 'Histone_modifier'
        gene_list.append('EP300')
        gene_types['JMJD3'] = 'Histone_modifier'
        gene_list.append('JMJD3')
        gene_types['BCL2'] = 'Survival'
        gene_list.append('BCL2')
        gene_types['BCL6'] = 'Survival'
        gene_list.append('BCL6')
        gene_types['BCL10'] = 'Survival'
        gene_list.append('BCL10')
        gene_types['TNFAIP3'] = 'Survival'
        gene_list.append('TNFAIP3')
        gene_types['REL'] = 'Survival'
        gene_list.append('REL')
        gene_types['FAS'] = 'Survival'
        gene_list.append('FAS')
        gene_types['MYC'] = 'Survival'
        gene_list.append('MYC')
        gene_types['SOCS1'] = 'Survival'
        gene_list.append('SOCS1')
        gene_types['CARD11'] = 'Survival'
        gene_list.append('CARD11')
        gene_types['MYD88'] = 'Survival'
        gene_list.append('MYD88')
        gene_types['CD97B'] = 'Survival'
        gene_list.append('CD79B')
        gene_types['PRDM1'] = 'Differentiation'
        gene_list.append('PRDM1')
        gene_types['IKZF3'] = 'Differentiation'
        gene_list.append('IKZF3')
        gene_types['IRF4'] = 'Differentiation'
        gene_types['IRF8'] = 'Differentiation'
        gene_list.append('IRF8')
        gene_types['HIST1H1C'] = 'Histone_core'
        gene_list.append('HIST1H1C')
        gene_types['HIST1H1E'] = 'Histone_core'
        gene_list.append('HIST1H1E')
        gene_types['HIST1H1B'] = 'Histone_core'
        gene_list.append('HIST1H1B')
        gene_types['HIST1H2AC'] = 'Histone_core'
        gene_list.append('HIST1H2AC')
        gene_types['HIST1H2AG'] = 'Histone_core'
        gene_list.append('HIST1H2AG')
        gene_types['HIST1H2AK'] = 'Histone_core'
        gene_types['HIST1H2BC'] = 'Histone_core'
        gene_types['HIST1H2BD'] = 'Histone_core'
        gene_list.append('HIST1H2BD')
        gene_types['HIST1H2BF'] = 'Histone_core'
        gene_types['HIST1H2BM'] = 'Histone_core'
        gene_types['HIST1H2BO'] = 'Histone_core'
        gene_types['HIST1H3A'] = 'Histone_core'
        gene_types['HIST1H3C'] = 'Histone_core'
        gene_types['HIST1H3D'] = 'Histone_core'
        gene_list.append('HIST1H3D')
        gene_types['HIST1H3F'] = 'Histone_core'
        gene_list.append('HIST1H3F')
        gene_types['HIST1H3H'] = 'Histone_core'
        gene_types['HIST1H4I'] = 'Histone_core'
        gene_types['HIST1H4K'] = 'Histone_core'
        gene_types['HIST2H2BE'] = 'Histone_core'
        gene_list.append('HIST2H2BE')
        gene_types['HIST2H2BF'] = 'Histone_core'
        gene_list.append('HIST2H2BF')
        gene_types['CCND3'] = 'Cell_cycle'
        gene_types['RB1'] = 'Cell_cycle'
        gene_list.append('RB1')
        gene_types['XPO1'] = 'mrna_shuttling'
        gene_types['ZFP36L1'] = 'mrna_shuttling'
        gene_types['TNPO1'] = 'mrna_shuttling'
        gene_types['TP53'] = 'Genome_instability'
        gene_list.append('TP53')
        gene_types['CDKN2A'] = 'Genome_instability'
        gene_list.append('CDKN2A')
        gene_types['CDKN2B'] = 'Genome_instability'
        gene_types['CD58'] = 'Immune'
        gene_list.append('CD58')
        gene_types['B2M'] = 'Immune'
        gene_list.append('B2M')
        gene_types['MIRN155'] = 'Other'
        gene_list.append('MIRN155')
        gene_types['MIRN17'] = 'Other'
        gene_list.append('MIRN17')
        ##gene_types['FOXP1'] = 'Other'
        gene_list.append('FOXP1')
        gene_types['GNA13'] = 'Other'
        gene_list.append('GNA13')
        gene_types['SGK1'] = 'Other'
        gene_list.append('SGK1')
        #gene_types['FOXO1'] = 'Other'
        #gene_list.append('FOXO1')
        #gene_types['PCLO'] = 'Other'
        #gene_list.append('PCLO')
        gene_types['TBL1XR1'] = 'Other'
        gene_types['TP63'] = 'Other'
        return (gene_list, gene_types)

class Sample():
    """Stores all useful details and functions for a sample"""
    def __init__(self, db_object, sample_id = None, sample_name = None):
        self.db = db_object
        cursor = db_object.cursor()
        if sample_name:
            query = "select sample.id, sample_id from sample where sample_id  = '%s'" % sample_name
        elif sample_id:
            query = "select sample.id, sample_id from sample where sample.id  = %s" % sample_id
        cursor.execute(query)
        result = cursor.fetchone()
        self.id, self.sample_name = result
    def __str__(self):
        return "%s %s" % (self.id,self.sample_name)
class Patient():
    """Stores all useful details and functions for a patient"""

    def __init__(self, db_object, patient_id = None, res_id = None, sample_name = None, library_name = None):
        self.db = db_object
        cursor = db_object.cursor()
        query = 'select sample.id, patient.id, res_id, sample.sample_id, sex from patient, sample, library where library.sample_id  = sample.id and patient.id = sample.patient_id and'
        if patient_id:
            query = query + ' patient.id = %s' % patient_id
        elif res_id:
            query = query + " res_id = '%s'"
        elif sample_name:
            query = query + " sample.sample_id = '%s'" % sample_name
        elif library_name:
            query = query + " library_name = '%s'" % library_name
        else:
            print 'no unique identifier provided'
            exit()
        #print query
        cursor.execute(query)
        result = cursor.fetchone()
        self.sample_id, self.id, self.res_id, self.sample_name, self.sex = result
    def __str__(self):
        return "%s %s %s" % (self.id,self.sample_name,self.res_id)
    def getLibrary(self, sample_type, library_type):
        """Retrieves the library for this patient of a given sample type (tumour/normal) and library type (genome/RNA-seq"""
        cursor = self.db.cursor()
        query = "select library.id from library, sample, patient where library.sample_id = sample.id and sample.patient_id = patient.id and library_type = '%s' and sample_type = '%s' and patient.id = %s"
        if sample_type == 'tumour' or sample_type == 'normal':
            pass
        else:
            print 'error, sample type %s not allowed' % sample_type
            exit()
        if library_type == 'genome' or library_type == 'exome' or library_type == 'RNA-seq':
            pass
        else:
            print 'error, library_type %s not allowed' % library_type
            exit()
        query = query % (library_type, sample_type, self.id)
        cursor.execute(query)
        result = cursor.fetchone()
        if result == None:
            return None
        library_id = result[0]
        library_obj = Library(self.db, library_id=library_id)
        return library_obj


class Library(cancerGenomeDB):
    """stores all useful details and functions for a library"""

    def __init__(self, db_object, library_id = None, library_name = None, sample_name = None, force_create = None):
        self.db = db_object
        #query = 'select library.id, library_name, library_type, sample.sample_id, sample.subtype, patient.res_id, patient.id, bam_location, reference_genome_fasta, average_coverage from sample, library, patient where patient.id = sample.patient_id and library.sample_id = sample.id and '
        query = 'select library.id, library_name, library_type, sample.sample_id, sample.subtype, patient.res_id, patient.id from sample, library, patient where patient.id = sample.patient_id and library.sample_id = sample.id and '
        if library_id:
            query = query + 'library.id = %i' % library_id
        elif library_name:
            query = query + "library_name = '%s'" % library_name
        else:
            print 'no library_id or library_name specified'
            exit()
        cursor = self.db.cursor()
        cursor.execute(query)
        data = cursor.fetchone()
        #print data
        if data:

            #self.id, library_name, ltype, sample_name, subtype, patient_external_id, patient_id, bam_location, ref_genome, average_coverage = data
            self.id, library_name, ltype, sample_name, subtype, patient_external_id, patient_id  = data
            self.library_name = library_name
            self.subtype = subtype
            self.library_type = ltype
            self.sample_name = sample_name
            self.patient_external_id = patient_external_id
            self.patient_id = patient_id
            #self.bam_location = bam_location
            #self.reference_genome = ref_genome
            #self.average_coverage = average_coverage

        elif(force_create):
            if not sample_name:
                print "error, need to supply sample information"
                exit()
            print "creating new library using %s and %s" % (library_name,sample_name)
            #first add sample record, if necessary
            query = "select id from sample where sample_id = '%s'" % sample_name
            cursor.execute(query)
            res = cursor.fetchone()
            try:
                sample_id = res[0]
            except IndexError:
                to_add = 1
            if to_add:
                query = "insert into sample (sample_id) values('%s')" % sample_name
                print query
        else:
            self.id = None
            self.library_name = None
            self.sample_name = None
            self.patient_id = None
    def __str__(self):
        return "%s %s" % (self.library_name,self.sample_name)
    def getSpliceSiteSNVs(self):
        return cancerGenomeDB._getSpliceSiteSNVs(self,library_id=self.id)
    def getIndels(self,limits=None):
        """Call the getIndels function specifying this library"""
        return cancerGenomeDB._getIndels(self, library_id=self.id,limits=limits)
    def countFusions(self):
        '''Just count the number of fusions in this library'''
        query = "select count(*) from fusion_transcript, event where library_id = %s and event.id = fusion_transcript.event_id" % self.id
        #print query
        cursor = self.db.cursor()
        cursor.execute(query)
        num = cursor.fetchone()[0]
        return(num)
    def getMutations(self, status = None, limits=None):
        """Call the getMutations function specifying this library"""
        if not limits:
            limits = {}
        if status:
            limits['status'] = status
        return cancerGenomeDB._getMutations(self, limits=limits, library_id=self.id)

    def setBamPath(self, bam_path):
        """Set the location of the most current bam file for this library, normally you don't have to call this"""
        cursor = self.db.cursor()
        query = "update library set bam_location = '%s' where id = %s" % (bam_path, self.id)
        print query
        cursor.execute(query)
    def updateCoverage(self):
        '''Set "average_coverage" in library table. Only needs to be done once, or whenever library bam file is updated with new data'''
        bam_path = self.bam_location
        stats_file = bam_path + "stats"
        stats = open(stats_file,'r')
        num = 0
        for line in stats:
            if line.startswith("Estimate_for_X_coverage"):
                line = line.rstrip("\n")
                (blah,num) = line.split()
                print num
        stats.close()
        stats = open(stats_file,'r')
        mapped = 0
        dupe = 0
        length = 0
        if num == 0:
            for line in stats:
                line = line.rstrip()
                if line.startswith("Read_length"):
                    (blah,length) = line.split()
                elif line.startswith("Number_of_Duplicates"):
                    (blah,dupe) = line.split()
                elif line.startswith("Number_Reads_Aligned"):
                    (blah,mapped) = line.split()
        if mapped > 0:
            num = (100 * (int(mapped) - int(dupe)))/2850000000
        if num > 0:
            query = "update library set average_coverage = %s where id = %s" %(num,self.id)
            cursor = self.db.cursor()
            cursor.execute(query)
    def __cmp__(self, other):
        if self.id == other.id:
            return 0
        elif self.id < other.id:
            return -1
        else:
            return 1
    def __hash__(self):
        return hash(self.id)

class Gene(cancerGenomeDB):
    """stores all useful details and functions for a gene"""

    def __init__(self, db_object, gene_id = None, ensembl_id = None, gene_symbol = None, gene_name= None):
        """collect useful details from the database for a gene based on its id or ensembl_id"""
        self.db = db_object
        query = 'select id, ensembl_id, gene_symbol, chromosome, start_position, end_position, biotype from gene where '
        if gene_id:
            query = query + 'id = %i' % gene_id
        elif ensembl_id:
            query = query + "ensembl_id = '%s'" % ensembl_id
        elif gene_symbol:
            query = query + "gene_symbol = '%s'" % gene_symbol
        elif gene_name:
            query = query + "gene_name = '%s'" % gene_name
        else:
            print 'No gene_id or ensembl_id specified'
            exit()
        cursor = self.db.cursor()
        #print query
        #exit()
        cursor.execute(query)
        data = cursor.fetchone()
        try:
            gene_id, ensembl_id, gene_symbol, chr, start, end, biotype = data
        except TypeError as e:
            print '%s not in database' % gene_symbol
            #print e
            #print query
            self.ensembl_id = None
            return

        self.id = gene_id
        self.gene_id = gene_id
        self.ensembl_id = ensembl_id
        self.gene_symbol = gene_symbol
        self.chromosome = chr
        self.chromosome_short = self.chromosome[3:]
        self.start = start
        self.end = end
        self.biotype = biotype
        self.distance_to_event = -1
    def __str__(self):
        return "%s %s %s %s %s %s" % (self.id,self.chromosome,self.start,self.end,self.ensembl_id,self.gene_symbol)
    def __cmp__(self, other):
        if self.id == other.id:
            return 0
        elif self.id < other.id:
            return -1
        else:
            return 1

    def __hash__(self):
        return hash(self.id)
    def getCDSLength(self):
        '''Calculate the length of the CDS for longest isoform'''
        transcripts = self.getTranscripts(longest=True)
        if len(transcripts) > 0:
            cds_length = transcripts[0].cds_length
        else:
            cds_length = 0
        return cds_length


    #def getSNVs(self, limits = None):
    #    limits['gene.id'] = '="%s"' % self.id
    #    return self._getMutations(limits=limits)
    def getExpression(self,sample_names=None,res_ids=None):
        '''Get the RPKM values for the gene in all (or some) libraries'''
        cursor = self.db.cursor()
        query = "select library_id, gene_expression.average_coverage, normalized_coverage from gene_expression, library, sample, patient where patient.id = sample.patient_id and gene_id = %s and sample.id = library.sample_id and library.id = gene_expression.library_id" % self.id
        if sample_names:
            query = query + " and sample.sample_id in ("
            for sample_name in sample_names:
                query = query + "'%s'," % sample_name
            query = query.rstrip(",")
            query = query + ")"
        elif res_ids:
            query = query + " and res_id in ("
            for res in res_ids:
                query = query + "'%s'," % res
            query = query.rstrip(",")
            query = query + ")"
        #print query
        cursor.execute(query)
        libs = []
        depth = []
        rpkm = []
        for row in cursor.fetchall():
            #print row
            lib = Library(self.db,library_id=row[0])
            libs.append(lib)
            depth.append(row[1])
            rpkm.append(row[2])
        return(libs,depth,rpkm)
    def getTranscripts(self, longest = None):
        """Create transcript objects for all the transcripts of this gene, return as list"""
        objects = []
        query = 'select transcript.id, cds_end - cds_start as cds_length from transcript where gene_id = %s order by cds_length DESC' % self.id
        cursor = self.db.cursor()
        cursor.execute(query)
        transcripts = cursor.fetchall()
        for transcript in transcripts:
            trans_obj = Transcript(self.db, transcript_id=transcript[0])
            objects.append(trans_obj)
            if longest:
                break

        return objects
    def getRegion(self):
        """Create a region object from the position of this gene"""
        region = Region(self.db,self.chromosome,self.start,self.end)
        return region
    def getAllSNVs(self,limits=None):
        '''Get all the SNVs in this gene from any library'''
        cursor = self.db.cursor()
        query = "select mutation.id from mutation, gene where gene.ensembl_id = mutation.gene and gene.id = %s" % self.id
        if limits:
            for key in limits:
                query = query + ' and %s %s' % (key, limits[key])
        cursor.execute(query)
        snvs = []
        try:
            for row in cursor.fetchall():
                snv_id = row[0]
                snv = SNV(self.db,snv_id)
                snvs.append(snv)
        except IndexError:
            print "nothing retrieved for gene %s" % self
            exit()
        return snvs
class Transcript(cancerGenomeDB):
    """Stores useful details for a transcript and has some functions for working with transcripts"""

    def __init__(self, db_object, transcript_id = None, ensembl_transcript_id = None):
        """initialize general details of object from database (transcript table)"""
        self.db = db_object
        cursor = self.db.cursor()
        if transcript_id:
            query = 'select transcript.ensembl_id, transcript.length, cds_start, cds_end, transcript.refseq, protein.id, gene.ensembl_id from gene, protein, transcript where gene.id = transcript.gene_id and protein.transcript_id = transcript.id and transcript.id = %s' % transcript_id
            cursor.execute(query)
            print query
            details = cursor.fetchone()
            print details
            try:
                self.id = transcript_id
            except TypeError:
                print "warning, no information for %s in database" % ensembl_transcript_id
                return None
            self.ensembl_id = details[0]

            self.length = details[1]
            self.cds_start = details[2]
            self.cds_end = details[3]
            self.cds_length = self.cds_end - self.cds_start + 1
            self.refseq = details[4]
            self.protein_id = details[5]
            self.ensembl_gene_id = details[6]
            self.gene = Gene(self.db, ensembl_id=self.ensembl_gene_id)
        else:
            query = "select transcript.id, transcript.length, cds_start, cds_end, transcript.refseq, gene.ensembl_id from transcript, gene where gene.id = transcript.gene_id and transcript.ensembl_id = '%s'" % ensembl_transcript_id
            print query
            cursor.execute(query)
            details = cursor.fetchone()
            print details
            try:
                self.id = details[0]
            except TypeError:
                print "warning, no information for %s in database" % ensembl_transcript_id
                return None
            self.ensembl_id = ensembl_transcript_id
            self.length = details[1]
            self.cds_start = details[2]
            self.cds_end = details[3]
            self.cds_length = self.cds_end - self.cds_start + 1
            self.refseq = details[4]
            self.ensembl_gene_id = details[5]
            self.gene = Gene(self.db, ensembl_id=self.ensembl_gene_id)


    def __str__(self):
        return "%s\t%s\t%s" % (self.ensembl_id,self.gene.gene_symbol,self.cds_length)
    def getSpliceSiteSNVs(self):
        '''Get the splice site SNVs in this gene that correspond to exons in this transcript'''
        all_in_gene = self._getSpliceSiteSNVs(gene=self.gene)
        all_in_transcript = []
        splice_sites_mutated = []
        exons = self.getExons()
        all_splice_site_positions = []
        for exon in exons:
            genome_start = exon.genome_start
            genome_end = exon.genome_end
            if exon.phase == -1:
                #UTR, ignore first splice site
                pass
            else:
                all_splice_site_positions.append(genome_start-1)
                all_splice_site_positions.append(genome_start-2)
            if exon.end_phase == -1:
                pass
            else:
                all_splice_site_positions.append(genome_end+1)
                all_splice_site_positions.append(genome_end+2)
        for splice_snv in all_in_gene:
            if splice_snv.position in all_splice_site_positions:
                all_in_transcript.append(splice_snv)
        return all_in_transcript

    def getGenomicPositions(self,coding_only=None):
        '''Return a list of all the genomic positions this transcript (or its CDS) correspond to'''
        exons = self.getExons()
        starts = []
        ends = []
        positions = []
        for exon in exons:
            if coding_only:
                starts.append(self.coding_genome_start)
                ends.append(self.coding_genome_end)
            else:
                starts.append(self.genome_start)
                ends.append(self.genome_end)
        for i in range(len(starts)):
            for pos in range(starts[i],ends[i]+1):
                positions.append(pos)
        return positions
    def getExons(self):
        """Create exon objects for all the transcripts of this gene, return as list"""
        objects = []
        query = 'select exon.id from exon where transcript_id = %s order by transcript_start' % self.id
        #print query
        cursor = self.db.cursor()
        cursor.execute(query)
        for result in cursor.fetchall():
            exon_ob = Exon(self.db, exon_id=result[0])
            if not hasattr(exon_ob,"chromosome_short"):
                continue
            objects.append(exon_ob)

        return objects
    def genomic2CDS(self,genomic_position,get_sequence=None,mutation=None):
        '''Convert a genomic position to a position relative to the CDS start of this transcript (accounting for introns etc) and also return the codon and codon position this site corresponds to'''
        exons = self.getExons()
        for exon in exons:
            #print "exon and %s" % genomic_position
            #print exon
            gstart = exon.genome_start
            gend = exon.genome_end
            if not (genomic_position >= gstart and genomic_position <= gend):
                print "skipping: %s-%s" % (gstart,gend)
                continue #wrong exon
            #print "MATCHING EXON"
            tstart = exon.transcript_start
            cds_start = self.cds_start + 1
            exon_cds_pos = tstart  - cds_start + 1#this is where in the CDS this exon starts
            print "strand: %s" % exon.strand
            if exon.strand == 1:
                offset = genomic_position - exon.genome_start+1
            else:
                offset = exon.genome_end - genomic_position+1
            if offset < 0:
                print "error, offset of %s, snv in UTR?" % offset
            this_cds_pos = offset + exon_cds_pos
            #depending on the position in the CDS, this will either be in position 1, 2 or 3 of the codon
            if this_cds_pos % 3:
                codon_pos = this_cds_pos % 3
            else:
                codon_pos = 3
            #now the tricky part, get the codon itself

            seq_offset = this_cds_pos - codon_pos
            compl = { "T": "A", "A": "T", "G": "C", "C": "G" }
            if get_sequence:
                seq = self.getSequence(cds_only=1)
                codon_seq = seq[seq_offset:seq_offset+3]
                large_seq = seq[seq_offset-42:seq_offset+45]
                from Bio.Seq import Seq
                seq_obj = Seq(large_seq)
                large_pep = str(seq_obj.translate())
                if mutation:
                    #replace base with mutant base, complement base if gene is on minus strand
                    if not exon.strand == 1:
                        mutation = compl[mutation]
                    print seq
                    print seq_offset
                    print codon_pos
                    print codon_seq
                    large_seq1 = seq[this_cds_pos-45-codon_pos:this_cds_pos-1]
                    large_seq2 = seq[this_cds_pos:this_cds_pos+45]
                    #large_seq = seq[seq_offset-42:seq_offset+45]
                    #this next bit of code will need to be modified to handle indels properly
                    if large_seq1 == "":
                        #error due to not enough bases, get the full sequence up to this position
                        large_seq1 = seq[:this_cds_pos-1]
                    if large_seq2 == "":
                        large_seq2 = seq[this_cds_pos:]
                    large_seq = large_seq1 +  seq[this_cds_pos-1:this_cds_pos] +  large_seq2
                    large_seq_mut = large_seq1 + mutation + large_seq2
                    print "mutant: %s\nnormal: %s" % (large_seq_mut,large_seq)
                    seq_obj = Seq(large_seq)
                    large_pep = str(seq_obj.translate())
                    mut_ob = Seq(large_seq_mut)
                    large_pep_mut = str(mut_ob.translate())
                    print "mutant: %s\nnormal: %s" % (large_pep_mut,large_pep)

                    return(large_seq,large_seq_mut,large_pep,large_pep_mut)
                else:
                    print "cds_pos: %s, codon: %s pos: %s (%s-%s) from %s\t%s and %s" % (this_cds_pos,codon_seq,codon_pos,seq_offset,seq_offset+3,genomic_position,large_seq,large_pep)
                    return(codon_pos,codon_seq,large_seq,large_pep)
            else:
                return(codon_pos,codon_seq)



        print "error, none of these exons seems to contain this SNV"
        print self
        return(-1,"")
    def cdna2genomic(self,cdna_position):
        '''Convert a cdna position to a coordinate position in the genome'''
        exons = self.getExons()
        for exon in exons:
            tstart = exon.transcript_start #where in the full-length transcript this exon begins
            tend = exon.transcript_end #where in the full-length transcript this exon ends
            if tstart <= cdna_position and tend >= cdna_position:
                pass
            else:
                continue
            gstart = exon.genome_start
            gend = exon.genome_end
            print exon
            if exon.strand == 1:
                exon_position = cdna_position - tstart
                genome_pos = gstart + exon_position
                print "%s from %s and %s" % (genome_pos,exon_position,cdna_position)
            else:
                exon_position = cdna_position - tstart + 1
                genome_pos = gend - exon_position + 1
                print exon
                print "%s from %s and %s" % (genome_pos,exon_position,cdna_position)
            return(genome_pos)
        return(-1,"")
    def getSNVs(self,limits=None):
        '''Get the SNVs that are within the CDS of this transcript'''
        cursor = self.db.cursor()
        #this code relied on codon table being populated for all isoforms, which doesn't appear to be the case
        #query = "select position, codon, codon_position from codon_gene where gene_id = %s " % self.gene.id
        #query = query + "and transcripts like '%" + self.ensembl_id + "%'"
        #print query
        #positions = {}
        #cursor.execute(query)
        #for res in cursor.fetchall():
        #    positions[res[0]] = (res[1],res[2])
        query = "select mutation.id from mutation, gene where gene.ensembl_id = mutation.gene and gene.id = %s" % self.gene.id
        if limits:
            for key in limits:
                query = query + ' and %s %s' % (key, limits[key])
        #print query
        cursor.execute(query)
        snvs = []
        for mutation_id in cursor.fetchall():
            mutation_obj = SNV(self.db,mutation_id)
            print mutation_obj
            mutation_pos = mutation_obj.position
            appended = 0
            (mutation_obj.codon_position,mutation_obj.codon) = self.genomic2CDS(mutation_pos)
            if mutation_obj.codon_position == -1:
                continue #not in the CDS of this isoform
            snvs.append(mutation_obj)
            appended = 1
        return snvs
    def getProtein(self):
        protein_obj = Protein(self.db, protein_id=self.protein_id)
        return protein_obj
    def getStrand(self):
        exons = self.getExons()
        strand = exons[0].strand
        return strand
    def getSequence(self,cds_only = None,genome_fasta=None):
        '''Using the exons, genome sequence and strand information, get the mRNA sequence'''
        genome_fasta='/projects/rmorin/common/genomes/hg19/GRCh37-lite.fa'
        exons = self.getExons()
        if len(exons) == 0:
            return ""
        mrna = ""
        for exon in exons:
            if not hasattr(exon,"gene"):
                continue
            if not hasattr(exon,"chromosome_short"):
                continue
            #print exon
            seq = exon.getSequence(coding=cds_only,genome_fasta=genome_fasta)  #will get either the CDS portion of exons or full length mRNA with UTR
            #print seq
            mrna = mrna + seq
        return mrna
    def calcTheoreticalGreenmanValues(self,genome_fasta=None):
        '''Using the sequence of the transcript and its codons, count how many of each mutation possibilities it has'''
        trans_seq = self.getSequence(cds_only=True,genome_fasta=genome_fasta)
        print trans_seq
        from Bio.Seq import Seq
        import re
        regex = re.compile("\w\w\w")
        codons = regex.findall(trans_seq)
        #print codons
        #exit()
        bases = ['A','T','G','C']
        l_ca = 0 # l == silent
        l_ct = 0
        l_cg = 0
        l_ta = 0
        l_tc = 0
        l_tg = 0
        m_ca = 0
        m_ct = 0
        m_cg = 0
        m_ta = 0
        m_tc = 0
        m_tg = 0
        n_ca = 0
        n_ct = 0
        n_cg = 0
        n_ta = 0
        n_tc = 0
        n_tg = 0
        for codon in codons:
            codon_seq = Seq(codon)
            aa = str(codon_seq.translate())

            for codon_pos in range(1,4):
                for other_base in bases:
                    base = codon[codon_pos-1]
                    if base == other_base:
                        continue
                    codon_new = list(codon)

                    codon_new[codon_pos-1] = other_base
                    codon_new_string = "".join(codon_new)
                    codon_seq_new = Seq(codon_new_string)
                    aa_new = str(codon_seq_new.translate())
                    if aa_new == aa:
                        #silent, L
                        #print "silent: %s to %s" % (aa,aa_new)
                        if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                            l_ca +=1
                        elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                            l_ct +=1
                        elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                            l_cg +=1
                        elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                            l_ta +=1
                        elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                            l_tc +=1
                        elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                            l_tg +=1
                    elif aa_new == "*" or (aa == "*" and not aa_new == "*"):
                        #truncating, N
                        #print "truncating: %s to %s" % (aa,aa_new)
                        if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                            n_ca +=1
                        elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                            n_ct +=1
                        elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                            n_cg +=1
                        elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                            n_ta +=1
                        elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                            n_tc +=1
                            #print "odd case: old: %s new: %s %s>%s (%s %s)" %(aa,aa_new,base,other_base,codon,codon_new)
                        elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                            n_tg +=1
                    else:
                        #missense
                        #print "missense: %s to %s" % (aa,aa_new)
                        if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                            m_ca +=1
                        elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                            m_ct +=1
                        elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                            m_cg +=1
                        elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                            m_ta +=1
                        elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                            m_tc +=1
                        elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                            m_tg +=1
        #ensembl L_ca L_cg L_ct L_ta L_tc L_tg M_ca M_cg M_ct M_ta M_tc M_tg N_ca N_cg N_ct N_ta N_tc N_tg S_ca S_cg S_ct S_ta S_tc S_tg, theoretical values (LMN)
        #also need S values
        splice_sites = []
        exons = self.getExons()
        if len(exons) > 0:
            import pysam
            reffile = pysam.Fastafile(genome_fasta)
            for exon in exons:
                print "exon"
                if exon.phase == -1:
                    #UTR, ignore first splice site
                    pass
                else:
                    seq = reffile.fetch(exon.chromosome_short, exon.genome_start - 3, exon.genome_start-1)
                    bases = list(seq)
                    for base in bases:
                        if base == "A":
                            base = "T"
                        elif base == "G":
                            base = "C"
                        splice_sites.append(base)
                if exon.end_phase == -1:
                    pass
                else:
                    seq = reffile.fetch(exon.chromosome_short, exon.genome_end , exon.genome_end+2)
                    bases = list(seq)
                    print seq
                    for base in bases:
                        if base == "A":
                            base = "T"
                        elif base == "G":
                            base = "C"
                        splice_sites.append(base)
                        print "base: %s" % base
        #count number of Cs (S_ca = S_ct = S_cg = #C)
        n_c_splice = 0
        n_t_splice = 0
        for base in splice_sites:
            if base == "T":
                n_t_splice +=1
            elif base == "C":
                n_c_splice +=1
            else:
                print "error on 1313"
                exit()
        return "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % (self.gene.ensembl_id,l_ca,l_cg,l_ct,l_ta,l_tc,l_tg,m_ca,m_cg,m_ct,m_ta,m_tc,m_tg,n_ca,n_cg,n_ct,n_ta,n_tc,n_tg,n_c_splice,n_t_splice)
    def calcObservedGreenmanValues(self,limits=None,from_file=None,coding_snv_file=None,splice_snv_file=None,gene_lookup=None):
        from Bio.Seq import Seq
        if from_file and gene_lookup:
            #for hg19, chandra's latest format
            gene_symbol_key = gene_lookup[self.gene.ensembl_id]
            #read contents of a file to calculate values
            print "HERE"

            handle = open(coding_snv_file,"r")
            #A1CF	10:52566514	+	Missense_Mutation	C/T/T	TCGA-46-3769-01A-01D-0983-08
            l_ca = 0 # l == silent
            l_ct = 0
            l_cg = 0
            l_ta = 0
            l_tc = 0
            l_tg = 0
            m_ca = 0
            m_ct = 0
            m_cg = 0
            m_ta = 0
            m_tc = 0
            m_tg = 0
            n_ca = 0
            n_ct = 0
            n_cg = 0
            n_ta = 0
            n_tc = 0
            n_tg = 0
            for line in handle:
                fields = line.split()
                if not fields[0] == gene_symbol_key:
                    continue
                mut_all = fields[4]
                vartype = fields[3]
                (base,other_base,other_other) = mut_all.split("/")
                print "change: %s to %s" % (base,other_base)
                other_base_copy = other_base
                #aa = fields[11]
                #aa_new = fields[12]
                strand = fields[2]
                if vartype == "Silent":
                    #silent, L
                    #print "silent: %s to %s" % (aa,aa_new)
                    if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                        l_ca +=1
                    elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                        l_ct +=1
                    elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                        l_cg +=1
                    elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                        l_ta +=1
                    elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                        l_tc +=1
                    elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                        l_tg +=1
                elif vartype == "Nonstop_Mutation" or vartype == "Nonsense_Mutation":
                    #truncating, N
                    #print "truncating: %s to %s" % (aa,aa_new)
                    if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                        n_ca +=1
                    elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                        n_ct +=1
                    elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                        n_cg +=1
                    elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                        n_ta +=1
                    elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                        n_tc +=1
                    elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                        n_tg +=1
                elif vartype == "Missense_Mutation":
                    #missense
                    print "missense"
                    if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                        m_ca +=1
                    elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                        m_ct +=1
                    elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                        m_cg +=1
                    elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                        m_ta +=1
                    elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                        m_tc +=1
                    elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                        m_tg +=1
            handle.close()
                #add any splice site mutations
            handle = open(splice_snv_file,'r')
            print "opened %s" % splice_snv_file
            S_c = 0
            S_t = 0
            #need to guess format (format 1)
            #    CDK5R1	ENSG00000176749	17:27839164	De_novo_Start_OutOfFrame	C	T	Somatic	exon2	c.C413T
            #   DYNC2H1	ENSG00000187240	11:102573871	Splice_Site	G	T	Somatic	.	.
            #(format 2)
            #chr11:64450957 G T G:2,T:2,0.0013232464,0.9986426453,0.0000341083,2 ENSG00000068971
            #new format:, hg19
            #A2M	12:9221440	+	Splice_Site	T/C/C	TCGA-39-5028-01A-01D-1441-08
            for line in handle:
                fields = line.split()
                vartype = fields[3]
                if vartype != "Splice_Site":
                    continue
                if not fields[0] == gene_symbol_key:
                    continue
                mut_all = fields[4]
                vartype = fields[3]
                (base,other_base,other_other) = mut_all.split("/")
                print "change: %s to %s" % (base,other_base)
                if base == "T" or base == "A":
                    S_t+=1
                elif base == "G" or base == "C":
                    S_c +=1

            handle.close()
            return "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % (self.gene.ensembl_id,l_ca,l_cg,l_ct,l_ta,l_tc,l_tg,m_ca,m_cg,m_ct,m_ta,m_tc,m_tg,n_ca,n_cg,n_ct,n_ta,n_tc,n_tg,S_c,S_t)
        elif from_file:
            #read contents of a file to calculate values
            handle = open(coding_snv_file,"r")
            #chr12:100826154 ENSG00000136048 1 GTC 3 134; ENST00000258534; C T null GTT V V SYNONYMOUS V134V
            l_ca = 0 # l == silent
            l_ct = 0
            l_cg = 0
            l_ta = 0
            l_tc = 0
            l_tg = 0
            m_ca = 0
            m_ct = 0
            m_cg = 0
            m_ta = 0
            m_tc = 0
            m_tg = 0
            n_ca = 0
            n_ct = 0
            n_cg = 0
            n_ta = 0
            n_tc = 0
            n_tg = 0
            for line in handle:
                fields = line.split()
                if len(fields) < 12:
                    continue
                if not fields[1] == self.gene.ensembl_id:
                    continue
                base = fields[7]
                other_base = fields[8]
                other_base_copy = other_base
                aa = fields[11]
                aa_new = fields[12]
                strand = fields[2]
                if aa_new == aa:
                    #silent, L
                    print "silent: %s to %s" % (aa,aa_new)
                    if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                        l_ca +=1
                    elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                        l_ct +=1
                    elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                        l_cg +=1
                    elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                        l_ta +=1
                    elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                        l_tc +=1
                    elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                        l_tg +=1
                elif aa_new == "*" or (aa == "*" and not aa_new == "*"):
                    #truncating, N
                    print "truncating: %s to %s" % (aa,aa_new)
                    if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                        n_ca +=1
                    elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                        n_ct +=1
                    elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                        n_cg +=1
                    elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                        n_ta +=1
                    elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                        n_tc +=1
                    elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                        n_tg +=1
                else:
                    #missense
                    print "missense: %s to %s" % (aa,aa_new)
                    if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                        m_ca +=1
                    elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                        m_ct +=1
                    elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                        m_cg +=1
                    elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                        m_ta +=1
                    elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                        m_tc +=1
                    elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                        m_tg +=1
            handle.close()
                #add any splice site mutations
            handle = open(splice_snv_file,'r')
            print "opened %s" % splice_snv_file
            S_c = 0
            S_t = 0
            #need to guess format (format 1)
            #    CDK5R1	ENSG00000176749	17:27839164	De_novo_Start_OutOfFrame	C	T	Somatic	exon2	c.C413T
            #   DYNC2H1	ENSG00000187240	11:102573871	Splice_Site	G	T	Somatic	.	.
            #(format 2)
            #chr11:64450957 G T G:2,T:2,0.0013232464,0.9986426453,0.0000341083,2 ENSG00000068971
            for line in handle:
                fields = line.split()
                #print fields
                if fields[0].startswith("chr"):
                    #format 2
                    if not fields[4] == self.gene.ensembl_id:
                        #print "skipping %s" % fields[4]
                        continue
                    base = fields[1]
                    other_base = fields[2]
                    print "splice affecting a %s" % base
                    if base == "T" or base == "A":
                        S_t+=1
                    elif base == "G" or base == "C":
                        S_c +=1
                else:
                    variant_type = fields[3]
                    #print "%s %s %s %s vs %s" % (variant_type, fields[7],fields[8],fields[1],self.gene.ensembl_id)
                    if not fields[1] == self.gene.ensembl_id:
                        continue
                    base = fields[4]
                    other_base = fields[5]
                    print "splice affecting a %s" % base
                    if base == "T" or base == "A":
                        S_t+=1
                    elif base == "G" or base == "C":
                        S_c +=1
            handle.close()
            return "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % (self.gene.ensembl_id,l_ca,l_cg,l_ct,l_ta,l_tc,l_tg,m_ca,m_cg,m_ct,m_ta,m_tc,m_tg,n_ca,n_cg,n_ct,n_ta,n_tc,n_tg,S_c,S_t)
        else:
            l_ca = 0 # l == silent
            l_ct = 0
            l_cg = 0
            l_ta = 0
            l_tc = 0
            l_tg = 0
            m_ca = 0
            m_ct = 0
            m_cg = 0
            m_ta = 0
            m_tc = 0
            m_tg = 0
            n_ca = 0
            n_ct = 0
            n_cg = 0
            n_ta = 0
            n_tc = 0
            n_tg = 0
            mutations = self.getSNVs(limits=limits)
            exons = self.getExons()
            strand = exons[0].strand
            for mut in mutations:
                base = mut.reference_base
                other_base = mut.nonreference_base
                other_base_copy = other_base
                if strand == -1:
                    if other_base == 'A':
                        other_base_copy = 'T'
                    elif other_base == 'C':
                        other_base_copy = 'G'
                    elif other_base == 'T':
                        other_base_copy = 'A'
                    elif other_base == 'G':
                        other_base_copy = 'C'

                codon = mut.codon
                codon_seq = Seq(codon)
                aa = str(codon_seq.translate())
                codon_pos = mut.codon_position
                codon_new = list(codon)
                codon_new[codon_pos-1] = other_base_copy
                codon_new_string = "".join(codon_new)
                codon_seq_new = Seq(codon_new_string)
                aa_new = str(codon_seq_new.translate())
                if aa_new == aa:
                    #silent, L
                    print "silent: %s to %s" % (aa,aa_new)
                    if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                        l_ca +=1
                    elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                        l_ct +=1
                    elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                        l_cg +=1
                    elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                        l_ta +=1
                    elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                        l_tc +=1
                    elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                        l_tg +=1
                elif aa_new == "*" or (aa == "*" and not aa_new == "*"):
                    #truncating, N
                    print "truncating: %s to %s" % (aa,aa_new)
                    if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                        n_ca +=1
                    elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                        n_ct +=1
                    elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                        n_cg +=1
                    elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                        n_ta +=1
                    elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                        n_tc +=1
                    elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                        n_tg +=1
                else:
                    #missense
                    print "missense: %s to %s" % (aa,aa_new)
                    if (base == 'C' and other_base == 'A') or (base == 'G' and other_base == 'T'):
                        m_ca +=1
                    elif (base == 'C' and other_base == 'T') or (base == 'G' and other_base == 'A'):
                        m_ct +=1
                    elif (base == 'C' and other_base == 'G') or (base == 'G' and other_base == 'C'):
                        m_cg +=1
                    elif (base == 'T' and other_base == 'A') or (base == 'A' and other_base == 'T'):
                        m_ta +=1
                    elif (base == 'T' and other_base == 'C') or (base == 'A' and other_base == 'G'):
                        m_tc +=1
                    elif (base == 'T' and other_base == 'G') or (base == 'A' and other_base == 'C'):
                        m_tg +=1
        #ensembl L_ca L_cg L_ct L_ta L_tc L_tg M_ca M_cg M_ct M_ta M_tc M_tg N_ca N_cg N_ct N_ta N_tc N_tg S_ca S_cg S_ct S_ta S_tc S_tg, theoretical values (LMN)
        #add any splice site mutations
        S_c = 0
        S_t = 0
        observed_splice_site_snvs = self.getSpliceSiteSNVs()
        for snv in observed_splice_site_snvs:
            ref = snv.base_change[:1]
            #print "%s reference is %s" % (snv.base_change,ref)
            #exit()
            if ref == "T" or ref == "A":
                S_t+=1
            elif ref == "G" or ref == "C":
                S_c +=1
        return "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % (self.gene.ensembl_id,l_ca,l_cg,l_ct,l_ta,l_tc,l_tg,m_ca,m_cg,m_ct,m_ta,m_tc,m_tg,n_ca,n_cg,n_ct,n_ta,n_tc,n_tg,S_c,S_t)
class Exon(cancerGenomeDB):
    """Stores useful details for an exon and has some functions for working with exons"""
    #genome_fasta='/projects/rmorin/common/genomes/hg19/GRCh37-lite.fa'
    def __init__(self, db_object, exon_id = None, ensembl_exon_id = None, transcript_id = None):
        """get general details of exon from database"""
        self.db = db_object
        if exon_id:
            query = 'select exon.id, ensembl_id, transcript_id, transcript_start, transcript_end, genome_start, genome_end, strand, phase, end_phase from exon where id = %s' % exon_id
        elif ensembl_exon_id and transcript_id:
            query = "select exon.id, ensembl_id, transcript_id, transcript_start, transcript_end, genome_start, genome_end, strand, phase, end_phase from exon where ensembl_id = '%s' and transcript_id = %s" % (ensembl_exon_id, transcript_id)
        else:
            print 'Error: need to supply both the ensembl_id and transcript_id'
            exit()
        cursor = self.db.cursor()
        cursor.execute(query)
        results = cursor.fetchone()
        self.id = results[0]
        self.ensembl_id = results[1]
        self.transcript_id = results[2]
        self.transcript_start = int(results[3])
        self.transcript_end = int(results[4])
        self.genome_start = int(results[5])
        self.genome_end = int(results[6])
        self.strand = int(results[7])
        self.phase = int(results[8])
        self.end_phase = int(results[9])
        #print query
        query = 'select gene.id, chromosome from gene, transcript where gene.id = transcript.gene_id and transcript.id = %s' % self.transcript_id
        #print query
        cursor.execute(query)
        self.gene_id, self.chromosome = cursor.fetchone()
        if self.chromosome == None:
            return None
        self.chromosome_short = self.chromosome[3:]
        self.gene = Gene(self.db, gene_id=self.gene_id)
        self.is_coding = 0
        if self.phase == -1 and self.end_phase == -1:
            self.is_coding = 1
        elif self.phase > -1:
            self.is_coding = 1
        elif self.end_phase > -1:
            self.is_coding = 1
        query = 'select cds_start, cds_end from transcript where id = %s' % self.transcript_id
        cursor.execute(query)
        result = cursor.fetchone()
        trans_cstart = result[0]
        trans_cend = result[1]
        if self.is_coding == 0:
            self.coding_start = -1
            self.coding_end = -1
        elif self.phase == -1:
            self.coding_start = trans_cstart - self.transcript_start + 1
            if self.end_phase == -1:
                #coding region ends in this exon
                self.coding_end = trans_cend - self.transcript_start + 1
            else:
                self.coding_end = self.transcript_end - self.transcript_start + 1
            #needs testing, probably has a bug
            if self.strand == 1:
                self.coding_genome_start = self.genome_start + trans_cstart - self.transcript_start + 1
                self.coding_genome_end = self.genome_end
            else:
                self.coding_genome_start = self.genome_start
                self.coding_genome_end = self.genome_end - (trans_cstart - self.transcript_start + 1)
        elif self.end_phase == -1:
            self.coding_start = 1
            self.coding_end = trans_cend - self.transcript_start + 1
            if self.strand == 1:
                self.coding_genome_start = self.genome_start
                self.coding_genome_end = self.genome_end - (trans_cend - self.transcript_start) + 1
            else:
                self.coding_genome_end = self.genome_end
                self.coding_genome_start = self.genome_start + (trans_cend - self.transcript_start) + 1
        else:
            self.coding_start = 1
            self.coding_end = self.transcript_end - self.transcript_start + 1
            self.coding_genome_start = self.genome_start
            self.coding_genome_end = self.genome_end
        if self.is_coding:
            self.coding_length = self.coding_end - self.coding_start + 1
            self.protein_length = float(self.coding_length) / 3
        else:
            self.coding_length = 0
            self.protein_length = 0

    def __cmp__(self, other):
        if self.id == other.id:
            return 0
        elif self.id < other.id:
            return -1
        else:
            return 1

    def __hash__(self):
        return hash(self.id)

    def __str__(self):
        try:
            to_return = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (self.id,
            self.gene.gene_symbol,
            self.chromosome,
            self.transcript_start,
            self.transcript_end,
            self.genome_start,
            self.genome_end,
            self.phase,
            self.end_phase,
            self.strand)
        except AttributeError as e:
            print e
            return "error here"
        return to_return

    def getSequence(self, coding = None,genome_fasta=None):
        """Using the strand, genome position, genome fasta file and phase information, get the sequence (if coding=yes, just the sequence corresponding to the cds)"""
        cursor = self.db.cursor()
        import pysam
        if self.chromosome == None:
            return ""
        print "using %s fasta" % genome_fasta
        #exit()
        reffile = pysam.Fastafile(genome_fasta)
        #print 'fetching: %s %s %s' % (self.chromosome_short, self.genome_start - 1, self.genome_end)
        full_seq = reffile.fetch(self.chromosome_short, self.genome_start - 1, self.genome_end)
        reffile.close()
        end_phase = self.end_phase
        phase = self.phase
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna, generic_protein
        seq_obj = Seq(full_seq, generic_dna)
        if self.strand == -1:
            #print 'revcomp!'
            seq_obj = seq_obj.reverse_complement()

        if coding:
            if phase == -1 and end_phase == -1: #single exon gene
                as_string = str(seq_obj)
                #print "%s to %s" % (self.coding_start,self.coding_end)
                return as_string[self.coding_start - 1:self.coding_end]
            elif phase == -1:
                as_string = str(seq_obj)
                return as_string[self.coding_start - 1:]
            elif end_phase == -1:
                as_string = str(seq_obj)
                return as_string[:self.coding_end]
            else:
                return str(seq_obj)
        else:
            return str(seq_obj)


class ProteinFeature():
    """stores useful details for an annotation/protein feature"""

    def __init__(self, db_object, feature_id):
        cursor = db_object.cursor()
        query = 'select protein_start, protein_end, source, description from protein_region where id = %s' % feature_id
        cursor.execute(query)
        start, end, source, description = cursor.fetchone()
        self.start = int(start)
        self.end = int(end)
        self.source = source
        self.description = description

    def __str__(self):
        to_return = '%s\t%s\t%s\t%s' % (self.start,
         self.end,
         self.source,
         self.description)
        return to_return


class Protein(cancerGenomeDB):
    """Stores useful details for a protein and has protein-centric functions"""

    def __init__(self, db_object, protein_id = None, ensembl_protein_id = None, transcript_id = None):
        """Initialize general details of object from database (protein table)"""
        self.db = db_object
        cursor = self.db.cursor()
        if protein_id:
            query = 'select protein.id, protein.ensembl_id, transcript_id, protein.length, gene_id from protein, transcript where transcript.id = protein.transcript_id and protein.id = %s' % protein_id
        elif transcript_id:
            query = 'select protein.id, protein.ensembl_id, transcript_id, protein.length, gene_id from protein, transcript where transcript.id = protein.transcript_id and transcript_id = %s' % transcript_id
        else:
            query = "select protein.id, protein.ensembl_id, transcript_id, protein.length, gene_id from protein, transcript where transcript.id = protein.transcript_id and protein.ensembl_id = '%s'" % ensembl_protein_id
        #print query
        cursor.execute(query)
        details = cursor.fetchone()
        self.id = details[0]
        self.ensembl_id = details[1]
        self.transcript_id = details[2]
        self.length = details[3]
        self.gene_id = details[4]
        self.gene = Gene(self.db, gene_id=self.gene_id)
        self.transcript = Transcript(self.db, transcript_id=self.transcript_id)

    def getExonLengths(self):
        """return list of exon lengths (in order)"""
        exon_objects = self.transcript.getExons()
        exon_lengths = []
        for exon in exon_objects:
            if exon.protein_length > 0:
                exon_lengths.append(exon.protein_length)

        return exon_lengths

    def getFeatures(self, protein_start = None, protein_end = None):
        """return position and details for all features annotated for this gene (e.g. domains). Remove or truncate features as necessary if protein_start or protein_end were provided"""
        cursor = self.db.cursor()
        features = []
        query = 'select id from protein_region where protein_id = %s order by source' % self.id
        cursor.execute(query)
        #print self.transcript.gene.gene_symbol, 'hello'
        for feature_id in cursor.fetchall():
            feature = ProteinFeature(self.db, feature_id[0])
            if protein_start:
                if feature.start >= protein_start:
                    pass
                elif feature.end > protein_start:
                    feature.start = protein_start
                else:
                    continue
            if protein_end:
                if feature.end <= protein_end:
                    pass
                elif feature.start < protein_end:
                    feature.end = protein_end
                else:
                    continue
            features.append(feature)

        return features

    def getMutations(self, mutation_class, status = None, library_id = None):
        """return position and details for all mutations in this gene in the database (eventually make more versatile), default is to return all (including status=unknown), which are from RNA-seq"""
        all_mutations = []
        cursor = self.db.cursor()
        if mutation_class == 'SNV':
            query = "select id, chromosome, position, base_change, annotation, protein_altering, library_id from mutation where gene = '%s'" % self.transcript.ensembl_gene_id
            if library_id:
                query += 'and library_id = %s' % library_id
            if status:
                status_string = ''
                for a in status:
                    status_string += "'%s'," % a

                status_string = status_string.rstrip(',')
                query += ' and status in (%s)' % status_string
            query = query + ' order by position'
            cursor.execute(query)
            for result in cursor.fetchall():
                mutation_id, chr, pos, change, annotations, prot_alt, library_id = result
                self.library_id = library_id
                query = "select transcripts, transcript_positions from codon_gene where chromosome = '%s' and position = %s and gene_id = %s" % (chr, pos, self.gene_id)
                cursor.execute(query)
                transcripts, positions = cursor.fetchone()
                transcripts = transcripts.rstrip(',')
                positions = positions.rstrip(',')
                annotations = annotations.rstrip(',')
                annotations_list = annotations.split(';')
                positions_list = positions.split(';')
                transcript_ids = transcripts.split(';')
                correct_position = None
                exit_loop = 0
                for i in range(len(annotations_list)):
                    if exit_loop:
                        break
                    transcript_id = transcript_ids[i]
                    codon_pos = positions_list[i]
                    annotation = annotations_list[i]
                    redundant_trans = transcript_id.split(',')
                    for trans_name in redundant_trans:
                        if trans_name == self.transcript.ensembl_id:
                            exit_loop = 1
                            correct_position = codon_pos
                            break
                        if exit_loop == 0:
                            continue

                if correct_position:
                    mutation_obj = ProteinMutation(self.db, self.transcript, correct_position, mutation_id, 'SNV', full_annotation=annotation, library_id=library_id)
                    all_mutations.append(mutation_obj)

            return all_mutations
        if mutation_class == 'indel':
            query = 'select indel.id, start, end, library_id from indel, event, gene_event where event.id = gene_event.event_id and event.id = indel.event_id and gene_id = %s' % self.gene_id
            if library_id:
                query += 'and library_id = %s' % library_id
            cursor = self.db.cursor()
            cursor.execute(query)
            for indel_data in cursor.fetchall():
                new_query = 'select transcripts, transcript_positions from codon_gene where gene_id = %s and position = %s' % (self.gene_id, indel_data[1])
                indel_id = indel_data[0]
                library_id = indel_data[3]
                cursor.execute(new_query)
                try:
                    trans, trans_pos = cursor.fetchone()
                except TypeError:
                    print 'indel not in this transcript'
                    break

                trans = trans.rstrip(';')
                trans_pos = trans_pos.rstrip(';')
                transcript_ids = trans.split(';')
                positions_list = trans_pos.split(';')
                exit_loop = 0
                correct_position = None
                for i in range(len(transcript_ids)):
                    if exit_loop:
                        break
                    transcript_id = transcript_ids[i]
                    codon_pos = positions_list[i]
                    redundant_trans = transcript_id.split(',')
                    for trans_name in redundant_trans:
                        if trans_name == self.transcript.ensembl_id:
                            exit_loop = 1
                            correct_position = codon_pos
                            break
                        if exit_loop == 0:
                            print 'no match found'
                            continue

                if correct_position:
                    mutation_obj = ProteinMutation(self.db, self.transcript, correct_position, indel_id, 'indel', library_id=library_id)
                    all_mutations.append(mutation_obj)

            return all_mutations

    def drawProteinImage(self):
        """Create and run an R script that draws the protein (exon structure) with all mutations and domains"""
        exon_lengths = self.getExonLengths()
        exon_centers = []
        last_end = 0
        for i in range(len(exon_lengths)):
            if i == 0:
                mid = float(exon_lengths[i]) / 2
                last_end = exon_lengths[i]
            else:
                mid = float(last_end) + float(exon_lengths[i]) / 2
                last_end = last_end + exon_lengths[i]
            exon_centers.append(mid)

        n_exons = len(exon_centers)
        image_end = last_end + 0.05 * last_end
        print 'dimensions = matrix(nrow=%i,ncol=2)' % n_exons
        lengths = ''
        centers = ''
        for i in range(len(exon_lengths)):
            lengths = lengths + '%f,' % exon_lengths[i]
            centers = centers + '%f,' % exon_centers[i]

        lengths = lengths.rstrip(',')
        centers = centers.rstrip(',')
        print 'dimensions[,1] = c(%s)' % lengths
        print 'dimensions[,2] = c(rep(0.06,%i))' % n_exons
        print 'centers = c(%s)' % centers
        print 'heights = rep(1,%s)' % n_exons
        print "symbols(centers,heights,rectangles=dimensions,inches=FALSE,bg=c('lightgrey','grey'),xlim=c(0,%s))" % image_end
        import math
        mutations = self.getMutations('SNV')
        circle_positions = []
        circle_y = []
        circle_cols = []
        for j in range(len(mutations)):
            mut = mutations[j]
            circle_positions.append(mut.position)
            colour = None
            if mut.type == 'synonymous':
                colour = 'blue'
            elif mut.type == 'non-synonymous':
                colour = 'green'
            elif mut.type == 'truncating' or mut.type == 'read-through':
                colour = 'red'
            circle_cols.append(colour)
            height = 1.1
            if j > 0:
                for k in range(len(circle_y)):
                    if math.fabs(mut.position - circle_positions[k]) < 50:
                        new_height = circle_y[k] + 0.01
                        if new_height > height:
                            height = new_height

                circle_y.append(height)
            else:
                circle_y.append(1.1)

        num_circles = len(circle_positions)
        circle_positions_string = ''
        circle_y_string = ''
        colour_string = ''
        l = len(circle_positions)
        for i in range(l):
            circle_positions_string = circle_positions_string + '%s,' % circle_positions[i]
            colour_string = colour_string + "'%s'," % circle_cols[i]
            circle_y_string = circle_y_string + '%s,' % circle_y[i]

        circle_positions_string = circle_positions_string.rstrip(',')
        colour_string = colour_string.rstrip(',')
        circle_y_string = circle_y_string.rstrip(',')
        radius = 0.003 * last_end
        circle_radius_string = 'rep(%s,%i)' % (radius, l)
        print 'symbols(c(%s),c(%s),circles=c(%s),bg=c(%s),add=TRUE,inches=FALSE)' % (circle_positions_string,
         circle_y_string,
         circle_radius_string,
         colour_string)
        indels = self.getMutations('indel')
        indel_positions = []
        indel_symbols = []
        height = 1.05
        for indel in indels:
            print '%s %s %s' % (indel.type, indel.position, indel.length)
            indel_positions.append(indel.position)
            if indel.type == 'ins':
                indel_symbols.append(25)
            else:
                indel_symbols.append(24)

        x_string = ''
        pch_string = ''
        height_string = ''
        k = 0
        num_ind = len(indels)
        for pos in range(num_ind):
            x_string = x_string + '%s,' % indel_positions[pos]
            height_string = height_string + '%s,' % height
            pch_string = pch_string + '%s,' % indel_symbols[k]
            k += 1

        height_string = height_string.rstrip(',')
        x_string = x_string.rstrip(',')
        pch_string = pch_string.rstrip(',')
        print "points(c(%s),c(%s),pch=c(%s),bg='orange')" % (x_string, height_string, pch_string)
        features = self.getFeatures()
        tracks = {}
        for feature in features:
            if feature.description == '':
                continue
            if tracks.has_key(feature.description):
                tracks[feature.description].append((feature.start, feature.end))
            else:
                tracks[feature.description] = [(feature.start, feature.end)]

        track_y = 0.95
        colours = ('gold', 'gold1', 'gold2', 'gold3', 'goldenrod', 'goldenrod1', 'goldenrod2', 'goldenrod3', 'deepskyblue', 'dodgerblue1', 'dodgerblue2')
        colind = 0
        text_offset = last_end + last_end / 100
        for description in tracks:
            centres = []
            widths = []
            regions = tracks[description]
            for r in regions:
                start = r[0]
                end = r[1]
                width = end - start + 1
                mid = float(start) + float(width) / 2
                centres.append(mid)
                widths.append(width)

            numb = len(widths)
            print 'dimensions = matrix(nrow=%i,ncol=2)' % numb
            centre_string = ''
            width_string = ''
            for i in range(numb):
                centre_string = centre_string + '%s,' % centres[i]
                width_string = width_string + '%s,' % widths[i]

            centre_string = centre_string.rstrip(',')
            width_string = width_string.rstrip(',')
            print 'dimensions[,1] = c(%s)' % width_string
            print 'dimensions[,2] = rep(0.012,%s)' % numb
            if colind >= len(colours):
                colind = 0
            print "symbols(c(%s),rep(%s,%s),rectangles=dimensions,bg='%s',inches=FALSE,add=TRUE)" % (centre_string,
             track_y,
             numb,
             colours[colind])
            print "text(%s,%s,'%s',cex=0.8,adj=c(0,0.5))" % (text_offset, track_y, description)
            track_y -= 0.015
            colind += 1

class CNV(cancerGenomeDB):
    '''Stores basic details for a CNV and cnv-specific functions'''
    def __init__(self,db_object,cnv_id):
        self.db = db_object
        cursor = self.db.cursor()
        self.id = cnv_id
        query = "select chromosome, segment_start, segment_end, segment_state, size, library_id, sample.sample_id, sample.id, event_id from cnv, event, library, sample where sample.id = library.sample_id and library.id = event.library_id and event.id = cnv.event_id and cnv.id = %s" % cnv_id
        cursor.execute(query)
        (self.chromosome,self.start,self.end,state,self.length,self.library_id,self.sample_name,self.sample_id,self.event_id) = cursor.fetchone()
        self.state = int(state)
    def getAffectedGenes(self,return_type='gene_object'):
        #this method was too slow, use gene_event table join instead!!!
        #genes = self.getGenesInRegion(self.chromosome,self.start,self.end)
        cursor = self.db.cursor()
        query = "select gene_id, relevant_target from gene_event where event_id = %s" % self.event_id
        cursor.execute(query)
        genes = {}
        try:
            for result in cursor.fetchall():
                gene_id = result[0]
                relevant_target = result[1]
                if return_type == 'gene_object':
                    gene = Gene(self.db,gene_id=gene_id)
                    genes[gene] = relevant_target
                elif return_type == "gene_id":
                    genes[gene_id] = relevant_target
        except TypeError:
            pass
        return genes
    def annotateAffectedGenes(self):
        '''Should only need to be called once per CNV (results are stored in database and can be accessed with getAffectedGenes)'''
        #select gene_symbol from gene where chromosome = 'chr7' and end_position >= 52912 and start_position <= 141694121
        cursor = self.db.cursor()
        query = "select id from gene where chromosome = '%s' and end_position >= %s and start_position <= %s" % (self.chromosome,self.start,self.end)
        cursor.execute(query)
        all_results = cursor.fetchall()
        if not all_results:
            pass
        else:
            for result in all_results:
                gene_id = result[0]
                check_query = "select count(*) from gene_event where event_id = %s and gene_id = %s" % (self.event_id,gene_id)
                #print check_query
                cursor.execute(check_query)
                counted = int(cursor.fetchone()[0])
                if counted:
                    pass
                else:
                    update_query = "insert into gene_event (gene_id,event_id) values(%s,%s)" % (gene_id,self.event_id)
                    #print update_query
                    cursor.execute(update_query)
    def __str__(self):
        return "%s:%s-%s\t%s\t%s" % (self.chromosome,self.start,self.end,self.length,self.state)
class Indel(cancerGenomeDB):
    """Stores basic details for an indel and indel-specific functions"""
    def __init__(self, db_object, indel_id):
        self.db = db_object
        self.id = indel_id
        #print "creating object for %s" % indel_id
        cursor = self.db.cursor()
        query = 'select gene.id, gene.chromosome, start, end, length, event_type, alt, event.to_validate, event.library_id, sample.sample_id, indel.validation_outcome from indel, event, gene_event, gene, library, sample where sample.id = library.sample_id and library.id = event.library_id and gene.id = gene_event.gene_id and event.id = indel.event_id and gene_event.event_id = event.id and indel.id = %s' % indel_id
        cursor.execute(query)
        self.gene_id = None
        #print query
        gene_obj = []
        first = 1
        for results in cursor.fetchall():
            print results
            (gene_id, chromosome, start, end, length, event_type, alt, to_validate, library_id, sample_name, validation_outcome) = results
            if first:
                self.chromosome = chromosome
                self.start = start
                self.end = end
                self.length = length
                self.type = event_type
                self.alt = alt
                self.library_id = library_id
                self.library = Library(self.db,library_id=library_id)
                self.sample_name = sample_name
                if to_validate == 'yes':
                    self.to_validate = True
                else:
                    self.to_validate = None
                self.validation_outcome = validation_outcome
            self.gene_id = gene_id
            gene_ob = Gene(db_object, gene_id=gene_id)
            gene_obj.append(gene_ob)
        #print results
        self.genes = gene_obj
        #now handle any indels with no record in gene_event table
        if not self.gene_id:
            query = 'select chromosome, start, end, alt, event.to_validate, event.library_id, sample.sample_id, indel.validation_outcome from indel, event, library, sample where sample.id = library.sample_id and library.id = event.library_id and event.id = indel.event_id and indel.id = %s' % indel_id
            cursor.execute(query)
            results = cursor.fetchone()
            (self.chromosome, self.start, self.end, self.alt, self.to_validate, self.library_id, self.sample_name, self.validation_outcome) = results
            self.library = Library(self.db,library_id=self.library_id)
    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.library.library_name,self.sample_name,self.type,self.id,self.chromosome,self.start,self.end, self.validation_outcome)
    def associatedCNV(self,window_size = 5000):
        '''Return boundaries of CNV that most likely corresponds to this indel'''
        cursor = self.db.cursor()
        query = "select count(*) from cnv, event where event.id = cnv.event_id and library_id = %s and (segment_start > %s - %s and segment_start < %s + %s) and (segment_end > %s - %s and segment_end < %s + %s) and segment_state < 2" % (self.library_id,self.start,window_size,self.start,window_size,self.end,window_size,self.end,window_size)
        #print query
        cursor.execute(query)
        n = cursor.fetchone()[0]
        print n
        if n > 0:
            query = "select cnv.id from cnv, event where event.id = cnv.event_id and library_id = %s and (segment_start > %s - %s and segment_start < %s + %s) and (segment_end > %s - %s and segment_end < %s + %s) and segment_state < 2" % (self.library_id,self.start,window_size,self.start,window_size,self.end,window_size,self.end,window_size)
            cursor.execute(query)
            cnv_id = cursor.fetchone()[0]
            print "CNV id: %s" % cnv_id
            return cnv_id
        else:
            return None
    def getEffectOnGenes(self):
        '''Return a hash of gene objects affected by indel, value is effect'''
        cursor = self.db.cursor()
        genes = self.getGenesInRegion(self.chromosome,self.start,self.end)
        gene_effect = {}
        for gene in genes:
            gene_start = gene.start
            gene_end = gene.end
            if self.start > gene_start and self.end < gene_end:
                query = "select exon.id from exon, transcript, gene where gene.id = transcript.gene_id and exon.transcript_id = transcript.id and gene.id = %s and ((genome_start > %s and genome_end < %s) or (genome_start < %s and genome_end > %s) or (genome_start < %s and genome_end > %s) or (genome_start < %s and genome_end > %s))" % (gene.id,self.start,self.end,self.start,self.end,self.end,self.end,self.start,self.start)
                #print query
                cursor.execute(query)
                exons = []
                try:
                    for result in cursor.fetchone():
                        #print result
                        exon = Exon(self.db,exon_id=result)
                        exons.append(exon)
                except TypeError:
                    pass
                gene_effect[gene] = exons
            else:
                gene_effect[gene] = "entire"
        return gene_effect
    def updateValidatorStatus(self, status):
        """Update the status for this indel, default is false since any indels not in the validator output are false positives"""
        query = "update indel set validator_result = '%s' where id= %s" % (status, self.id)
        cursor = self.db.cursor()
        cursor.execute(query)
class Region():
    '''Generic region, for anything not in database that we want to know mutation status, cnv status etc'''
    def __init__(self,db_object,chromosome,start,end=None,name=None):
        self.db = db_object
        self.chromosome = chromosome
        self.start = start
        if end:
            self.end  = end
        else:
            self.end = start
        if name:
            self.name = name
    def __str__(self):
        string_version = "%s\t%s\t%s\t%s" % (self.name,self.chromosome,self.start,self.end)
        return string_version
    def getSomaticSNVs(self,library=None):
        '''get all the somatic SNVs within a given region (organized by individual library)'''
        #query based on position and then create SNV objects for each
        #db_object, mutation_id, genome_wide=True
        cursor = self.db.cursor()
        query = "select somatic_mutation.id from somatic_mutation, event where event.id = somatic_mutation.event_id and chromosome = '%s' and position >= %s and position <= %s" % (self.chromosome,self.start,self.end)

        if library:
            query = query + " and library_id = %s" % library.id
        #print query
        cursor.execute(query)
        snv_ids = []
        for result in cursor.fetchall():
            snv_ids.append(result[0])
        snvs = []
        for snv_id in snv_ids:
            snv = SNV(self.db,snv_id,genome_wide=True)
            snvs.append(snv)
        return snvs
    def getCopyNumberStatus(self,library):
        """load/return the state from the CNV table for this mutation"""
        cursor = self.db.cursor()
        query = "select segment_state, size from cnv, event where cnv.event_id = event.id and event.library_id = %s and chromosome = '%s' and ((segment_start <= %s and segment_end >= %s) or (segment_start >= %s and segment_start <= %s) or (segment_end >= %s and segment_end <= %s))" % (library.id,self.chromosome,self.start,self.end,self.start,self.end,self.start,self.end)
        #print query
        cursor.execute(query)
        cnv_types = []
        cnv_sizes = []
        try:
            state, size = cursor.fetchone()
        except TypeError:
            return ("NEUT",2)
        #assumes that if no segment is found, there is no event there
        if state == 1:
            cnv_type = 'HETD'
        elif state == 0:
            cnv_type = 'HOMD'
        elif state == 2:
            cnv_type = 'NEUT'
        elif state == 3:
            cnv_type = 'GAIN'
        elif state == 4:
            cnv_type = 'AMP'
        else:
            cnv_type = 'HAMP'
        self.cnv_state = state
        self.cnv_type = cnv_type
        return (cnv_type,state)


class SpliceSiteSNV():
    '''Stores basic details for a splice site SNV'''
    def __init__(self,db_object,splice_site_snv_id):
        self.db = db_object
        cursor = db_object.cursor()
        query = 'select gene_id, event.id, chromosome, position, base_change, validation_outcome, library_name, sample.sample_id from sample, gene_event, splice_site_snv, event, library where sample.id = library.sample_id and library.id = event.library_id and splice_site_snv.event_id = event.id and gene_event.event_id = event.id and splice_site_snv.id = %s' % splice_site_snv_id
        #print query
        cursor.execute(query)
        (gene_id,event_id,chrom,pos,base_change,validation_outcome,library_name,sample_name) = cursor.fetchone()
        gene_obj = Gene(db_object,gene_id=gene_id)
        self.gene = gene_obj
        self.gene_id = gene_id
        self.id = splice_site_snv_id
        self.event_id = event_id
        self.chromosome = chrom
        self.position = pos
        self.base_change = base_change
        self.validation_outcome = validation_outcome
        self.library_name = library_name
        self.sample_name= sample_name
        self.reference_base = base_change[0]
        self.nonreference_base = base_change[2]
    def __str__(self):
        return "%s %s %s %s %s" % (self.id, self.chromosome, self.position, self.base_change, self.library_name)
class SNV(cancerGenomeDB):
    """Stores basic details for a single base change"""
    def __init__(self, db_object, mutation_id, genome_wide=None):
        self.db = db_object
        cursor = db_object.cursor()
        if genome_wide:
            query = 'select somatic_mutation.id, library_id, somatic_mutation.chromosome, position, base_change, validation_outcome, sample.sample_id from event, somatic_mutation, library, sample  where sample.id = library.sample_id and library.id = event.library_id and somatic_mutation.event_id = event.id and somatic_mutation.id = %s' % mutation_id
            cursor.execute(query)
            self.id, library_id, self.chromosome, self.position, base_change, self.validation_outcome, self.sample_name = cursor.fetchone()
            protein_altering = None
            self.annotation = None
        else:
            query = 'select mutation.gene, gene.id, mutation.chromosome, position, annotation, base_change, protein_altering, library_id, sample.sample_id, status, validation_outcome, validation_method from mutation, gene, library, sample where sample.id = library.sample_id and library.id = mutation.library_id and gene.ensembl_id = mutation.gene and mutation.id = %s' % mutation_id
            cursor.execute(query)
            self.ensembl_id, self.gene_id, self.chromosome, self.position, self.annotation, base_change, protein_altering, library_id, self.sample_name, self.status, self.validation_outcome, self.validation_method = cursor.fetchone()
            self.id = mutation_id
        self.library = Library(self.db, library_id=library_id)
        self.chromosome_short = self.chromosome[3:]
        self.base_change = base_change
        if protein_altering == 'yes':
            self.type = 'non-synonymous'
        elif protein_altering == None:
            self.type = 'non-genic'
        else:
            self.type = 'synonymous'
        self.reference_base = base_change[0]
        self.nonreference_base = base_change[2]
    def __cmp__(self, other):
        if self.id == other.id:
            return 0
        elif self.id < other.id:
            return -1
        else:
            return 1
    def asMAF(self,for_mutsig = None,fake_id=None):
        '''return MAF formatted version of variant'''
        maf_fields = {}
        #dummy values and defaults
        if fake_id:
            maf_fields['Entrez_Gene_Id'] = fake_id
        else:
            maf_fields['Entrez_Gene_Id'] = 1
        maf_fields['Center'] = "bcgsc.ca"
        maf_fields['NCBI_Build'] = 'hg18'
        maf_fields['BAM_file'] = 'fake'
        maf_fields['Strand'] = "+"
        maf_fields['Tumor_Seq_Allele2'] = ""
        maf_fields['Genome_Change'] = self.base_change
        maf_fields['dbSNP_RS'] = 'novel'
        maf_fields['dbSNP_Val_Status']=''
        maf_fields['Variant_Type'] = 'SNP'
        if hasattr(self,'gene_id'):
            gene_obj = Gene(self.db, gene_id=self.gene_id)
            if gene_obj.gene_symbol == "":
                gene_obj.gene_symbol = gene_obj.gene_id
            maf_fields['Hugo_Symbol'] = gene_obj.gene_symbol
        else:
            maf_fields['Hugo_Symbol'] = "-"
        #print maf_fields['Hugo_Symbol']
        maf_fields['Chromosome'] = self.chromosome
        maf_fields['Start_position'] = self.position
        maf_fields['End_position'] = self.position
        annos = self.annotation.split(";")
        anno1 = annos[0]
        if self.type== 'non-synonymous':
            maf_fields['is_coding'] = 1
            maf_fields['is_silent'] = 0
            if anno1.endswith("*"):
                maf_fields['Variant_Classification'] = 'Nonsense_Mutation'
            else:
                maf_fields['Variant_Classification'] = 'Missense_Mutation'
        elif self.type == 'synonymous':
            maf_fields['is_silent'] = 1
            maf_fields['is_coding'] = 1
            maf_fields['Variant_Classification'] = 'Silent'
        else:
            maf_fields['is_coding'] = 0
            maf_fields['is_silent'] = 1
            maf_fields['Variant_Classification'] = "IGR"
        (ref,mut) = self.base_change.split(">")
        maf_fields['Annotation_Transcript'] = "none"
        maf_fields['Transcript_Strand'] = "+"
        maf_fields['Transcript_Exon'] = 1
        maf_fields['Transcript_Position'] = 1
        maf_fields['cDNA_Change'] = 1
        maf_fields['Codon_Change'] = 1
        maf_fields['Protein_Change'] = self.annotation
        maf_fields['Reference_Allele'] = ref
        maf_fields['Tumor_Seq_Allele1'] = mut
        maf_fields['Tumor_Sample_Barcode'] = self.sample_name + "-Tumor"
        maf_fields['Matched_Norm_Sample_Barcode'] = self.sample_name + "-Normal"
        maf_fields['Match_Norm_Seq_Allele1'] = ""
        maf_fields['Match_Norm_Seq_Allele2'] = ""
        maf_fields['Tumor_Validation_Allele1'] = ''
        maf_fields['Tumor_Validation_Allele2'] = ''
        maf_fields['Match_Norm_Validation_Allele1'] = ''
        maf_fields['Match_Norm_Validation_Allele2'] = ''
        maf_fields['Verification_Status']= ''
        maf_fields['Validation_Status']= ''
        maf_fields['Mutation_Status']= ''
        maf_fields['Sequencing_Phase']= ''
        maf_fields['Sequence_Source']= ''
        maf_fields['Validation_Method'] = ""
        maf_fields['Score']= ''
        maf_fields['BAM_File']= ''
        maf_fields['Sequencer']= ''
        maf_fields['Tumor_Sample_UUID']= ''
        maf_fields['Matched_Norm_Sample_UUID']= ''
        #ordered = ['Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome','Start_Position','End_Position','Strand','Variant_Classification','Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','dbSNP_RS','dbSNP_Val_Status','Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode','Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1','Tumor_Validation_Allele2','Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status','Validation_Status','Mutation_Status','Sequencing_Phase','Sequence_Source','Validation_Method','Score','BAM_File','Sequencer','Tumor_Sample_UUID','Matched_Norm_Sample_UUID']
        ordered = ['Hugo_Symbol',
                   'Entrez_Gene_Id',
                   'Center',
                   'NCBI_Build',
                   'Chromosome',
                   'Start_position',
                    'End_position',
                    'Strand',
                    'Variant_Classification',
                    'Variant_Type',
                    'Reference_Allele',
                    'Tumor_Seq_Allele1',
                    'Tumor_Seq_Allele2',
                    'dbSNP_RS',
                    'dbSNP_Val_Status',
                    'Tumor_Sample_Barcode',
                    'Matched_Norm_Sample_Barcode',
                    'Match_Norm_Seq_Allele1',
                    'Match_Norm_Seq_Allele2',
                    'Tumor_Validation_Allele1',
                    'Tumor_Validation_Allele2',
                    'Match_Norm_Validation_Allele1',
                    'Match_Norm_Validation_Allele2',
                    'Verification_Status',
                    'Validation_Status',
                    'Mutation_Status',
                    'Sequencing_Phase',
                    'Sequence_Source',
                    'Validation_Method',
                    'Score',
                    'BAM_file',
                    'Sequencer',
                    'Genome_Change',
                    'Annotation_Transcript',
                    'Transcript_Strand',
                    'Transcript_Exon',
                    'Transcript_Position',
                    'cDNA_Change',
                    'Codon_Change',
                    'Protein_Change']
        if for_mutsig:
            ordered.append("is_coding")
            ordered.append("is_silent")
            ordered.append("categ")
        #need to calculate the category (1. CpG transitions, 2. CpG transversions, 3. C:G transitions (CG>TA), 4. C:G transversions, 5. A:T transitions, 6. A:T transversions, 7. null+indel mutations) for each
        #get CpG context for C and G refbase
        category = 0
        import pysam
        genome_fasta = '/projects/rmorin/common/genomes/human_all.fasta'
        reffile = pysam.Fastafile(genome_fasta)
        #seq_right = reffile.fetch(self.chromosome_short, self.position,self.position+1)
        #seq_left = reffile.fetch(self.chromosome_short,self.position-1,self.position)
        triplet = reffile.fetch(self.chromosome_short,self.position-2,self.position+1)
        print triplet
        (seq_left,check_ref,seq_right) = list(triplet)
        if not ref == check_ref:
            print "error comparing %s to %s" % (triplet,ref)
        if maf_fields['Variant_Classification'] == 'Nonsense_Mutation':
            category = 7
        elif ref == "C":
            print "left: %s self: %s right: %s " % (seq_left,ref,seq_right)
            if seq_right == "G":
                print "CpG for %s" % self.position
                if mut == "T":
                    #transition
                    category = 1
                else:
                    category = 2
            else:
                print "not CpG: %s %s" % (self.position,ref)
                if mut == "T":
                    #transition
                    category = 3
                else:
                    category = 4

        elif ref == "G":
            print "left: %s self: %s right: %s " % (seq_left,ref,seq_right)
            if seq_left == "C":
                print "CpG for %s" % self.position
                if mut == "A":
                    #transition
                    category = 1
                else:
                    category = 2
            else:
                print "not CpG: %s %s" % (self.position,ref)
                if mut == "A":
                    #transition
                    category = 3
                else:
                    category = 4
        elif ref == "A" or ref == "T":
            if mut == "G" and ref == "A":
                category = 5
            elif ref == "A":
                category = 6
            elif ref == "T" and mut == "C":
                category = 5
            elif ref == "T":
                category = 6
        else:
            category = 7
        print "category: %s" % category
        maf_fields['categ'] = category

        maf_return = ""
        maf_header = ""
        printed = 0
        for key in ordered:
            if not maf_fields.has_key(key):
                print "key missing: %s"  % key
                exit()
            val = maf_fields[key]
            if not printed:
                maf_header = str(key)
                maf_return = str(val)
                printed = 1
            else:
                maf_return = maf_return + "\t" + str(val)
                maf_header = maf_header + "\t" + key
        #print self
        #print maf_header
        #print maf_return
        return(maf_header,maf_return)
    def __hash__(self):
        return hash(self.id)
    def __str__(self):
        return '%s:%s,%s,%s,%s,%s %s' % (
            self.chromosome,
            self.position,
            self.base_change,
            self.annotation,
            self.ensembl_id,
            self.validation_outcome,self.library)
    def getAffectedTranscripts(self):
        '''Return a list of transcript objects that have this SNV in one of their exons'''
        query = "select distinct transcript.id from exon, gene, transcript, mutation where gene = '%s' and mutation.gene = gene.ensembl_id and gene.id = transcript.gene_id and exon.transcript_id = transcript.id and exon.genome_start>= %s and exon.genome_end >= %s" % (self.ensembl_id,self.position,self.position)
        cursor = self.db.cursor()
        cursor.execute(query)
        trans_ob = []
        for row in cursor.fetchall():
            trans_obj = Transcript(self.db,transcript_id=row[0])
            print trans_obj
            trans_ob.append(trans_obj)
        return trans_ob
    def getCopyNumberStatus(self):
        #I think this doesn't work yet. Need the query to include the position of the mutaiton!
        """load/return the state from the CNV table for this mutation"""
        cursor = self.db.cursor()
        query = 'select segment_state, size from cnv, event where cnv.event_id = event.id and event.library_id = %s' % self.library.id
        cursor.execute(query)
        state, size = cursor.fetchone()
        if state == 1:
            cnv_type = 'HETD'
        elif state == 0:
            cnv_type == 'HOMD'
        elif state == 2:
            cnv_type = 'NEUT'
        elif state == 3:
            cnv_type = 'GAIN'
        elif state == 4:
            cnv_type = 'AMP'
        else:
            cnv_type = 'HAMP'
        self.cnv_state = state
        self.cnv_type = cnv_type
        return cnv_type
    def getLohValues(self):
        cursor = self.db.cursor()
        #check whether library is RNA-seq, if so, need to get corresponding genome library
        lib_type = self.library.library_type
        if lib_type == "RNA-seq":
            glib_res = self.getLibraries(type='genome',sample_name=self.library.sample_name,sample_type='tumour')
            try:
                glib = glib_res[0]
                print "fetched genome lib %s for RNA-seq library %s" % (self.library,glib)
            except KeyError:
                print "no result in glib_res"
                exit(0)
            query = "select loh.id, copy_number, loh_state from loh, event where event.id = loh.event_id and event.type = 'loh' and library_id = %s and chromosome = '%s' and segment_start <= %s and segment_end >= %s" % (glib.id,self.chromosome,self.position,self.position)
        else:
            query = "select loh.id, copy_number, loh_state from loh, event where event.id = loh.event_id and event.type = 'loh' and library_id = %s and chromosome = '%s' and segment_start <= %s and segment_end >= %s" % (self.library.id,self.chromosome,self.position,self.position)
        print query
        cursor.execute(query)
        try:
            (loh_id,cn,state) = cursor.fetchone()
        except TypeError:
            print "no loh segment for %s" % self.id
            return(-1,-1)
        print "for %s, id is %s, cn is %s and state is %s" % (self.id,loh_id,cn,state)
        loh = 0
        if state == "DLOH" or state == "NLOH" or state == "ALOH":
            loh = 1
        return(cn,loh)
    def getReseqValues(self):
        cursor = self.db.cursor()
        #reseq_high_quality_reference, reseq_high_quality_nonreference
        query = "select reseq_high_quality_reference, reseq_high_quality_nonreference from mutation where id = %s" % self.id
        cursor.execute(query)
        result = cursor.fetchone()
        return(result)
    def updateValitionDetails(self, field, new_value):
        """Allow updating of details concering the mutations validation status (or to_validate status)"""
        cursor = self.db.cursor()
        query = "update mutation set %s = '%s' where id = %s"
        if field == 'to_validate':
            if new_value == 'yes' or new_value == 'no':
                pass
            else:
                print 'value %s not allowed for %s' % (new_value, field)
                return 0
        elif field == 'status':
            if new_value == 'unknown' or new_value == 'valid_somatic' or new_value == 'valid_germline' or new_value == 'artefact':
                pass
            else:
                print 'value %s not allowed for %s' % (new_value, field)
                return 0
        elif field == 'validation_outcome':
            if new_value == 'somatic' or new_value == 'germline' or new_value == 'false' or new_value == 'unclear' or new_value == 'known_somatic':
                pass
            else:
                print 'value %s not allowed for %s' % (new_value, field)
                return 0
        elif field == 'validation_method':
            if new_value == 'sanger' or new_value == 'illumina' or new_value == 'iontorrent':
                pass
            else:
                print 'value %s not allowed for %s' % (new_value, field)
                return 0
        query = query % (field, new_value, self.id)
        cursor.execute(query)

    def updateBaseCounts(self, ref_count = None, nonref_count = None, total = None, rnaseq_ref_count = None, rnaseq_nonref_count = None, germline_nonreference_base_count = None, high_quality_germline_coverage = None, reseq_high_quality_reference = None, reseq_high_quality_nonreference = None,reseq_high_quality_reference_germline = None, reseq_high_quality_nonreference_germline = None):
        """Used to update reference_base_count and nonreference_base_count and high_quality_coverage in mutation table"""
        cursor = self.db.cursor()
        if isinstance(ref_count, int):
            query = 'update mutation set reference_base_count = %s where id = %s' % (ref_count, self.id)
            cursor.execute(query)
        if isinstance(nonref_count, int):
            query = 'update mutation set nonreference_base_count = %s where id = %s' % (nonref_count, self.id)
            cursor.execute(query)
        if isinstance(total, int):
            query = 'update mutation set high_quality_coverage = %s where id = %s' % (total, self.id)
            cursor.execute(query)
        if isinstance(rnaseq_ref_count, int):
            query = 'update mutation set rnaseq_reference_base_count = %s where id = %s' % (rnaseq_ref_count, self.id)
            cursor.execute(query)
        if isinstance(rnaseq_nonref_count, int):
            query = 'update mutation set rnaseq_nonreference_base_count = %s where id = %s' % (rnaseq_nonref_count, self.id)
            cursor.execute(query)
        if isinstance(high_quality_germline_coverage, int):
            query = 'update mutation set high_quality_germline_coverage = %s where id = %s' % (high_quality_germline_coverage, self.id)
            cursor.execute(query)
        if isinstance(germline_nonreference_base_count, int):
            query = 'update mutation set germline_nonreference_base_count = %s where id = %s' % (germline_nonreference_base_count, self.id)
            cursor.execute(query)
        if isinstance(reseq_high_quality_reference, int):
            query = 'update mutation set reseq_high_quality_reference = %s where id = %s' % (reseq_high_quality_reference, self.id)
            cursor.execute(query)
        if isinstance(reseq_high_quality_nonreference, int):
            query = 'update mutation set reseq_high_quality_nonreference = %s where id = %s' % (reseq_high_quality_nonreference, self.id)
            cursor.execute(query)
        if isinstance(reseq_high_quality_reference_germline, int):
            query = 'update mutation set reseq_high_quality_reference_germline = %s where id = %s' % (reseq_high_quality_reference_germline, self.id)
            cursor.execute(query)
        if isinstance(reseq_high_quality_nonreference_germline, int):
            query = 'update mutation set reseq_high_quality_nonreference_germline = %s where id = %s' % (reseq_high_quality_nonreference_germline, self.id)
            cursor.execute(query)

class ProteinMutation(cancerGenomeDB):
    """stores useful details for a mutation (relative to a given protein/transcript)"""

    def __init__(self, db_object, transcript, amino_acid_position, id, mutation_class, full_annotation = None, library_id = None):
        cursor = db_object.cursor()
        if mutation_class == 'SNV':
            query = 'select id from transcript_mutation where mutation_id = %s and transcript_id = %s' % (id, transcript.id)
            cursor.execute(query)
            try:
                test = cursor.fetchone()[0]
            except TypeError:
                print 'must add this to transcript_mutation table %s %s' % (id, transcript.id)
                query = "insert into transcript_mutation (mutation_id,transcript_id,annotation) values(%s,%s,'%s')" % (id, transcript.id, full_annotation)
                cursor.execute(query)

        if library_id:
            self.library = Library(db_object, library_id=library_id)
        if mutation_class == 'SNV':
            m = re.search('([A-Z*])([0-9]+)([A-Z*])', full_annotation)
            aa1 = m.group(1)
            blah = m.group(2)
            aa2 = m.group(3)
            type = None
            if aa1 == aa2:
                type = 'synonymous'
            elif aa2 == '*':
                type = 'truncating'
            elif aa1 == '*':
                type = 'read-through'
            else:
                type = 'non-synonymous'
            self.annotation = full_annotation
        elif mutation_class == 'indel':
            query = 'select start, end, length, event_type from indel where id = %s' % id
            cursor.execute(query)
            start, end, length, type = cursor.fetchone()
            self.length = length
        self.position = int(amino_acid_position)
        self.transcript = transcript
        self.id = id
        self.type = type


class FusionTranscript(cancerGenomeDB):
    """Allows hybrid transcripts containing regions from 2 transcripts to be handled and analyzed, annotated etc.
    Warning, this class is incomplete, buggy, and doesn't yet incude the contig sequences, which are important for it to actually derive the CDS properly
    Use at your own risk. You have been warned!"""
    def __init__(self, db, fusion_transcript_id,reverse_order=None):
        self.db = db
        cursor = self.db.cursor()
        query = "select description, genomic_regions, strands, mutation_type, library_id from event, fusion_transcript where event.id = fusion_transcript.event_id and fusion_transcript.id = %s" % fusion_transcript_id
        cursor.execute(query)
        (description,genomic_regions,strands,mutation_type,lib_id) = cursor.fetchone()
        #print strands
        #print description
        library = Library(self.db,library_id=lib_id)
        self.library = library
        (part1,part2) = genomic_regions.split(",")
        if reverse_order:
            part1, part2 = part2, part1
        (chromosome1,region1) = part1.split(":")
        (chromosome2,region2) = part2.split(":")
        try:
            start1, end1 = region1.split('-')
            start2, end2 = region2.split('-')
        except ValueError:
            print "this fusion is not formatted properly"
            return None
        self.start1 = int(start1) - 1
        self.start2 = int(start2) - 1
        self.end1 = int(end1)
        self.end2 = int(end2)
        try:
            (strand1,strand2) = strands.split(",")
        except ValueError:
            print "this fusion is not formatted properly"
            return None
        except AttributeError:
            print "something wrong here too"
            return None
        if reverse_order:
            strand1, strand2  = strand2, strand1
        #if strand1 == '-':
        #    end1 = self.start1 + 100
        #if strand2 == '-':
        #    start2 = self.end2 - 100

        genes1 = self.getGenesInRegion(chromosome1, self.start1, self.end1)
        n1 = len(genes1)
        genes2 = self.getGenesInRegion(chromosome2, self.start2, self.end2)
        n2 = len(genes2)
        #print '%s and %s genes' % (n1, n2)

        #print '%s\n%s' % (genes1[0].gene_symbol, genes2[0].gene_symbol)
        if len(genes1) > 0 and hasattr(genes1[0],"id"):
            gene1 = genes1[0]
            self.gene1 = gene1
        if len(genes2) > 0 and hasattr(genes2[0],"id"):
            gene2 = genes2[0]
            self.gene2 = gene2
        if n1 < 1 or n2 < 1:
            print "no genes found for %s %s %s %s %s %s" % (chromosome1,self.start1,self.end1,chromosome2,self.start2,self.end2)
            return None
        if not hasattr(self,"gene1") and hasattr(self,"gene2"):
            return None
        transcripts_1 = gene1.getTranscripts(longest=1)
        transcripts_2 = gene2.getTranscripts(longest=1)

        if len(transcripts_1) < 1 or len(transcripts_2)< 1:
            return None
        else:
            #print "both have transcripts"
            trans1 = transcripts_1[0]
            trans2 = transcripts_2[0]
            #print trans1
            #print "and"
            #print trans2

        ex1 = trans1.getExons()
        self.transcript1 = trans1
        self.transcript2 = trans2
        #print "transcript 1 exons:"
        #print ex1
        #print start1
        strand1 = ex1[0].strand
        fusion_exons = []
        for exon in ex1:
            start = exon.genome_start
            end = exon.genome_end
            if strand1 == -1:
                if start >= self.start1:
                    fusion_exons.append(exon)
            else:
                if start <= self.end1:
                    fusion_exons.append(exon)
        nex = len(fusion_exons)
        #print "%s exons from gene 1" % nex
        self.transcript1_last_exon = len(fusion_exons) - 1
        self.transcript2_first_exon = len(fusion_exons)
        ex2 = trans2.getExons()
        strand2 = ex2[0].strand
        #print "transcript 2 exons, strand is %s, gene is %s:" % (strand2,ex2[0].gene.gene_symbol)
        #print ex2
        #print start2
        added = 0
        for exon in ex2:
            start = exon.genome_start
            end = exon.genome_end
            if strand2 == -1:
                #print "minus strand"
                if start <= self.end2:
                    fusion_exons.append(exon)
                    added+=1
            else:
                if end >= self.end2:
                    fusion_exons.append(exon)
                    added+=1
        if added == 0:
            print "error, 0 exons added from gene2 using criteria strand: %s, end: %s" % (strand2,self.end2)
            #for exon in ex2:
                #print exon
            return None

        nex = len(fusion_exons)
        #print "%s exons after adding gene 2" % nex
        self.fusion_exons = fusion_exons
        full_seq = ''
        for ex in fusion_exons:
            #print ex
            #print '%s %s %s %s' % (ex.gene.gene_symbol,ex.genome_start,ex.genome_end,ex.id)
            coding_seq = ex.getSequence(coding=1)
            coding_start = ex.coding_start
            coding_end = ex.coding_end
            #print coding_seq, coding_start, coding_end
            full_seq += coding_seq

        self.exons = fusion_exons
        #print full_seq
        self.cds = full_seq
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna, generic_protein
        seq_obj = Seq(full_seq, generic_dna)
        peptide = seq_obj.translate()
        #print peptide
        self.protein1 = Protein(self.db, transcript_id=self.transcript1.id)
        self.protein2 = Protein(self.db, transcript_id=self.transcript2.id)
        self.transcript1_end_position = fusion_exons[self.transcript1_last_exon].transcript_end
        #print 'transcript1 end: %s and cds_start: %s' % (self.transcript1_end_position, self.transcript1.cds_start)
        #print fusion_exons[0]
        self.transcript2_start_position = fusion_exons[self.transcript2_first_exon].transcript_start
        #print self.transcript1_last_exon
        #print self.transcript2_first_exon
        #print fusion_exons[self.transcript1_last_exon]
        #print fusion_exons[self.transcript2_first_exon]
        #print 'transcript2 start: %s and cds_start2 = %s' % (self.transcript2_start_position, self.transcript2.cds_start)
        self.protein1_end_position = (self.transcript1_end_position - self.transcript1.cds_start) / 3 + 1
        self.protein2_start_position = (fusion_exons[self.transcript2_first_exon].transcript_start - self.transcript2.cds_start) / 3 + 1
        #print 'protein1 goes from 1-%s and protein2 goes from %s to end' % (self.protein1_end_position, self.protein2_start_position)
        features1 = self.protein1.getFeatures(protein_start=1, protein_end=self.protein1_end_position)
        all_features = []
        for feature in features1:
            all_features.append(feature)

        features2 = self.protein2.getFeatures(protein_start=self.protein2_start_position)
        for index, feature in enumerate(features2):
            feature.start -= self.protein2_start_position
            feature.end -= self.protein2_start_position
            feature.start += self.protein1_end_position
            feature.end += self.protein1_end_position
            all_features.append(feature)

        self.all_features = all_features

    def getExonLengths(self):
        """return list of exon lengths (in order)"""
        exon_objects = self.fusion_exons
        exon_lengths = []
        for exon in exon_objects:
            if exon.protein_length > 0:
                exon_lengths.append(exon.protein_length)

        return exon_lengths

    def drawProteinImage(self):
        """Create and run an R script that draws the protein (exon structure) with all mutations and domains"""
        exon_lengths = self.getExonLengths()
        print exon_lengths
        exon_centers = []
        last_end = 0
        for i in range(len(exon_lengths)):
            if i == 0:
                mid = float(exon_lengths[i]) / 2
                last_end = exon_lengths[i]
            else:
                mid = float(last_end) + float(exon_lengths[i]) / 2
                last_end = last_end + exon_lengths[i]
            exon_centers.append(mid)

        n_exons = len(exon_centers)
        image_end = last_end + 0.05 * last_end
        print 'dimensions = matrix(nrow=%i,ncol=2)' % n_exons
        lengths = ''
        centers = ''
        for i in range(len(exon_lengths)):
            lengths = lengths + '%f,' % exon_lengths[i]
            centers = centers + '%f,' % exon_centers[i]

        lengths = lengths.rstrip(',')
        centers = centers.rstrip(',')
        print 'dimensions[,1] = c(%s)' % lengths
        print 'dimensions[,2] = c(rep(0.06,%i))' % n_exons
        print 'centers = c(%s)' % centers
        print 'heights = rep(1,%s)' % n_exons
        print "symbols(centers,heights,rectangles=dimensions,inches=FALSE,bg=c('lightgrey','grey'),xlim=c(0,%s))" % image_end
        features = self.all_features
        tracks = {}
        for feature in features:
            if feature.description == '':
                continue
            if tracks.has_key(feature.description):
                tracks[feature.description].append((feature.start, feature.end))
            else:
                tracks[feature.description] = [(feature.start, feature.end)]

        track_y = 0.95
        colours = ('gold', 'gold1', 'gold2', 'gold3', 'goldenrod', 'goldenrod1', 'goldenrod2', 'goldenrod3', 'deepskyblue', 'dodgerblue1', 'dodgerblue2')
        colind = 0
        text_offset = last_end + last_end / 100
        for description in tracks:
            centres = []
            widths = []
            regions = tracks[description]
            for r in regions:
                start = r[0]
                end = r[1]
                width = end - start + 1
                mid = float(start) + float(width) / 2
                centres.append(mid)
                widths.append(width)

            numb = len(widths)
            print 'dimensions = matrix(nrow=%i,ncol=2)' % numb
            centre_string = ''
            width_string = ''
            for i in range(numb):
                centre_string = centre_string + '%s,' % centres[i]
                width_string = width_string + '%s,' % widths[i]

            centre_string = centre_string.rstrip(',')
            width_string = width_string.rstrip(',')
            print 'dimensions[,1] = c(%s)' % width_string
            print 'dimensions[,2] = rep(0.012,%s)' % numb
            if colind >= len(colours):
                colind = 0
            print "symbols(c(%s),rep(%s,%s),rectangles=dimensions,bg='%s',inches=FALSE,add=TRUE)" % (centre_string,
             track_y,
             numb,
             colours[colind])
            print "text(%s,%s,'%s',cex=0.8,adj=c(0,0.5))" % (text_offset, track_y, description)
            track_y -= 0.015
            colind += 1


class Rearrangement(cancerGenomeDB):
    """Stores useful details for a genomic rearrangement (which comprises two separate breakpoints)"""

    def __init__(self, db_object, event_id):
        cursor = db_object.cursor()
        query = 'select library_id, event.type, chromosome, position, manual_annotation, contigs, gene_id, sample.sample_id, genomic_break.validation_outcome from genomic_break, event, library, sample where sample.id = library.sample_id and library.id = event.library_id and event.id = event_id and event_id = %s order by position' % event_id
        cursor.execute(query)
        self.db = db_object
        chromosomes = []
        positions = []
        annotations = []
        for result in cursor.fetchall():
            library_id, event_type, chromosome, position, anno, contigs, gene_id, sample_name, validation_outcome = result
            chromosomes.append(chromosome)
            positions.append(int(position))
            annotations.append(anno)
            if hasattr(self, 'library'):
                pass
            else:
                self.library = Library(db_object, library_id=library_id)
                self.type = event_type

        self.annotations = annotations
        self.chromosomes = chromosomes
        self.sample_name = sample_name
        self.positions = positions
        self.event_id = event_id
        self.annotated_gene = gene_id
        self.type = event_type
        self.library_id = library_id
        self.contigs = contigs
        self.validation_outcome = validation_outcome
    def __str__(self):
        '''BED format for easy printing into a BED file for display in IGV, UCSC etc'''
        out_string1 = "%s\t%s\t%s\t%s\t%s\t%s" % (self.library.library_name,self.library.sample_name,self.type,self.contigs,self.chromosomes[0],self.positions[0])
        out_string2 = "%s\t%s\t%s\t%s\t%s" % (self.chromosomes[1],self.positions[1],self.event_id,self.validation_outcome,self.annotations[0])
        return out_string1 + "\t" + out_string2
    def associatedCNV(self,window_size = 5000):
        '''Return boundaries of CNV that most likely corresponds to this deletion'''
        if not self.type == "deletion":
            print "this only works on deletions"
            return None

        cursor = self.db.cursor()
        if self.positions[0] < self.positions[1]:
            start = self.positions[0]
            end = self.positions[1]
        else:
            start = self.positions[1]
            end = self.positions[0]
        query = "select count(*) from cnv, event where event.id = cnv.event_id and library_id = %s and (segment_start > %s - %s and segment_start < %s + %s) and (segment_end > %s - %s and segment_end < %s + %s) and segment_state < 2" % (self.library_id,start,window_size,start,window_size,end,window_size,end,window_size)
        #print query
        cursor.execute(query)
        n = cursor.fetchone()[0]
        print n
        if n > 0:
            query = "select cnv.id from cnv, event where event.id = cnv.event_id and library_id = %s and (segment_start > %s - %s and segment_start < %s + %s) and (segment_end > %s - %s and segment_end < %s + %s) and segment_state < 2" % (self.library_id,start,window_size,start,window_size,end,window_size,end,window_size)
            cursor.execute(query)
            cnv_id = cursor.fetchone()[0]
            print "CNV id: %s" % cnv_id
            return cnv_id
        else:
            return None
    def getContainedGenes(self, max_size=500000, symbol_required = None):
        """Get all genes within an event for duplications and deletions"""
        if self.type != 'deletion' and self.type != 'duplication':
            print 'Error, this function cannot be used'
            return []
        import math
        event_size = math.fabs(self.positions[0] - self.positions[1])
        if event_size > max_size:
            print 'Error, this event is larger than max_size %s' % event_size
            return []
        if self.positions[0] < self.positions[1]:
            start = self.positions[0]
            end = self.positions[1]
        else:
            start = self.positions[1]
            end = self.positions[0]
        all_genes = self.getGenesInRegion(self.chromosomes[0], start, end)
        some_genes = []
        if symbol_required:
            for gene in all_genes:
                if gene.gene_symbol:
                    some_genes.append(gene)

            return some_genes
        else:
            return all_genes

    def getProximalGenes(self, max_dist, max_genes, symbol_required = None):
        """Get determine all genes within max_dist of this rearrangement. Return a list of genes for each breakpoint, if more than max_genes exist, return nearest max_genes"""
        genes1 = self.getGenesInRegion(self.chromosomes[0], self.positions[0] - max_dist, self.positions[0] + max_dist)
        gene_dists = {}
        import math
        for gene in genes1:
            dist1 = math.fabs(gene.start - self.positions[0])
            dist2 = math.fabs(gene.end - self.positions[0])
            if dist1 > dist2:
                gene_dists[gene] = dist2
            else:
                gene_dists[gene] = dist1
            if gene.start < self.positions[0] and gene.end > self.positions[0]:
                gene_dists[gene] = 0

        genes_sorted_1 = []
        for key, value in sorted(gene_dists.iteritems(), key=lambda (k, v): (v, k)):
            if symbol_required:
                if key.gene_symbol == '':
                    continue
            gene_ob = key
            gene_ob.distance_to_event = int(value)
            genes_sorted_1.append(gene_ob)
            if len(genes_sorted_1) == max_genes:
                break

        genes2 = self.getGenesInRegion(self.chromosomes[1], self.positions[1] - max_dist, self.positions[1] + max_dist)
        gene_dists = {}
        for gene in genes2:
            dist1 = math.fabs(gene.start - self.positions[1])
            dist2 = math.fabs(gene.end - self.positions[1])
            if dist1 > dist2:
                gene_dists[gene] = dist2
            else:
                gene_dists[gene] = dist1
            if gene.start < self.positions[1] and gene.end > self.positions[1]:
                gene_dists[gene] = 0

        genes_sorted_2 = []
        for key, value in sorted(gene_dists.iteritems(), key=lambda (k, v): (v, k)):
            if symbol_required:
                if key.gene_symbol == '':
                    continue
            gene_ob = key
            gene_ob.distance_to_event = int(value)
            genes_sorted_2.append(gene_ob)
            if len(genes_sorted_2) == max_genes:
                break

        return (genes_sorted_1, genes_sorted_2)
