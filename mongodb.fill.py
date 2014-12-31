
from mongoalchemy.session import Session
from mongoalchemy.document import Document, Index 
from mongoalchemy.fields import *
import sqlite3
import json
DB = 'MTB_database_merged.db'

config_ref = {

	'reference': 
	{
		'snp_file':'dbsnp_137.hg19.dbsnp.only_3_first_cols.vcf',
		'gene_file':'genes.ref.csv',
		'DB':'../MTB_database_merged.db'
	},	
}
config_MTB = {
	'project_name': 'MTB',
	'genotypes': {
		'file': 'genotypes.MTB.csv',
		'fields': {
			'snp_id': 0,
			'snp_chr': 1,
			'snp_position': 2,
			'geno1': 3,
			'geno2': 4,
			'geno3': 5,
			'conditions': [
				'Stimulated',
				'Not Stimulated'				
			],
			'individuals': [
				{
					'name': 'Mtb2',
					'index': 6
				},
				{
					'name': 'Mtb3',
					'index': 7
				},
				{
					'name': 'Mtb5',
					'index': 8
				},
				{
					'name': 'Mtb6',
					'index': 9
				}
				,
				{
					'name': 'Mtb7',
					'index': 10
				}
				,
				{
					'name': 'Mtb8',
					'index': 11
				}
				,
				{
					'name': 'Mtb9',
					'index': 12
				}
				,
				{
					'name': 'Mtb10',
					'index': 13
				}
				,
				{
					'name': 'Mtb11',
					'index': 14
				}
				,
				{
					'name': 'Mtb12',
					'index': 15
				}
				,
				{
					'name': 'Mtb13',
					'index': 16
				}
				,
				{
					'name': 'Mtb14',
					'index': 17
				}
				,
				{
					'name': 'Mtb16',
					'index': 18
				}
				,
				{
					'name': 'Mtb17',
					'index': 19
				}
				,
				{
					'name': 'Mtb20',
					'index': 20
				}
				,
				{
					'name': 'Mtb22',
					'index': 21
				}
				,
				{
					'name': 'Mtb23',
					'index': 22
				}
				,
				{
					'name': 'Mtb24',
					'index': 23
				}
				,
				{
					'name': 'Mtb25',
					'index': 24
				}
				,
				{
					'name': 'Mtb26',
					'index': 25
				}
				,
				{
					'name': 'Mtb27',
					'index': 26
				}
				,
				{
					'name': 'Mtb28',
					'index': 27
				}
				,
				{
					'name': 'Mtb29',
					'index': 28
				}
				,
				{
					'name': 'Mtb30',
					'index': 29
				}
				,
				{
					'name': 'Mtb31',
					'index': 30
				}
				,
				{
					'name': 'Mtb32',
					'index': 31
				}
				,
				{
					'name': 'Mtb33',
					'index': 32
				}
				,
				{
					'name': 'Mtb38',
					'index': 33
				}
				,
				{
					'name': 'Mtb39',
					'index': 34
				}
				,
				{
					'name': 'Mtb40',
					'index': 35
				},
				{
					'name': 'Mtb41',
					'index': 36
				},
				{
					'name': 'Mtb42',
					'index': 37
				},
				{
					'name': 'Mtb43',
					'index': 38
				},
				{
					'name': 'Mtb44',
					'index': 39
				}
				,
				{
					'name': 'Mtb45',
					'index': 40
				}
				,
				{
					'name': 'Mtb46',
					'index': 41
				}
				,
				{
					'name': 'Mtb47',
					'index': 42
				}
				,
				{
					'name': 'Mtb49',
					'index': 43
				}
				,
				{
					'name': 'Mtb50',
					'index': 44
				}
				,
				{
					'name': 'Mtb51',
					'index': 45
				}
				,
				{
					'name': 'Mtb52',
					'index': 46
				}
				,
				{
					'name': 'Mtb53',
					'index': 47
				}
				,
				{
					'name': 'Mtb54',
					'index': 48
				}
				,
				{
					'name': 'Mtb55',
					'index': 49
				}
				,
				{
					'name': 'Mtb57',
					'index': 50
				}
				,
				{
					'name': 'Mtb59',
					'index': 51
				}
				,
				{
					'name': 'Mtb68',
					'index': 52
				}
				,
				{
					'name': 'Mtb69',
					'index': 53
				}
				,
				{
					'name': 'Mtb70',
					'index': 54
				}
				,
				{
					'name': 'Mtb71',
					'index': 55
				}
				,
				{
					'name': 'Mtb72',
					'index': 56
				}
				,
				{
					'name': 'Mtb74',
					'index': 57
				}
				,
				{
					'name': 'Mtb75',
					'index': 58
				}
				,
				{
					'name': 'Mtb77',
					'index': 59
				}
				,
				{
					'name': 'Mtb78',
					'index': 60
				}
				,
				{
					'name': 'Mtb80',
					'index': 61
				}
				,
				{
					'name': 'Mtb81',
					'index': 62
				}
				,
				{
					'name': 'Mtb82',
					'index': 63
				}
				,
				{
					'name': 'Mtb83',
					'index': 64
				}
				,
				{
					'name': 'Mtb84',
					'index': 65
				},
				{
					'name': 'Mtb85',
					'index': 66
				},
				{
					'name': 'Mtb86',
					'index': 67
				},
				{
					'name': 'Mtb87',
					'index': 68
				},
				{
					'name': 'Mtb88',
					'index': 69
				},
				{
					'name': 'Mtb90',
					'index': 70
				}

			]
		}
	},
	'expressions': {
		'file': [

		{
			'name':'expression.not_stim.MTB.csv',
			'condition':'Not Stimulated'
		},
		{
			'name':'expression.stim.MTB.csv',
			'condition':'Stimulated' 
		}
	],
		'fields': {
			'ensembl_gene_id': 0,
			'external_gene_id': 1,
			'individuals': [
				{
					'name': 'Mtb2',
					'index': 2
				},
				{
					'name': 'Mtb3',
					'index': 3
				},
				{
					'name': 'Mtb5',
					'index': 4
				},
				{
					'name': 'Mtb6',
					'index': 5
				}
				,
				{
					'name': 'Mtb7',
					'index': 6
				}
				,
				{
					'name': 'Mtb8',
					'index': 7
				}
				,
				{
					'name': 'Mtb9',
					'index': 8
				}
				,
				{
					'name': 'Mtb10',
					'index': 9
				}
				,
				{
					'name': 'Mtb11',
					'index': 10
				}
				,
				{
					'name': 'Mtb12',
					'index': 11
				}
				,
				{
					'name': 'Mtb13',
					'index': 12
				}
				,
				{
					'name': 'Mtb14',
					'index': 13
				}
				,
				{
					'name': 'Mtb16',
					'index': 14
				}
				,
				{
					'name': 'Mtb17',
					'index': 15
				}
				,
				{
					'name': 'Mtb20',
					'index': 16
				}
				,
				{
					'name': 'Mtb22',
					'index': 17
				}
				,
				{
					'name': 'Mtb23',
					'index': 18
				}
				,
				{
					'name': 'Mtb24',
					'index': 19
				}
				,
				{
					'name': 'Mtb25',
					'index': 20
				}
				,
				{
					'name': 'Mtb26',
					'index': 21
				}
				,
				{
					'name': 'Mtb27',
					'index': 22
				}
				,
				{
					'name': 'Mtb28',
					'index': 23
				}
				,
				{
					'name': 'Mtb29',
					'index': 24
				}
				,
				{
					'name': 'Mtb30',
					'index': 25
				}
				,
				{
					'name': 'Mtb31',
					'index': 26
				}
				,
				{
					'name': 'Mtb32',
					'index': 27
				}
				,
				{
					'name': 'Mtb33',
					'index': 28
				}
				,
				{
					'name': 'Mtb38',
					'index': 29
				}
				,
				{
					'name': 'Mtb39',
					'index': 30
				}
				,
				{
					'name': 'Mtb40',
					'index': 31
				},
				{
					'name': 'Mtb41',
					'index': 32
				},
				{
					'name': 'Mtb42',
					'index': 33
				},
				{
					'name': 'Mtb43',
					'index': 34
				},
				{
					'name': 'Mtb44',
					'index': 35
				}
				,
				{
					'name': 'Mtb45',
					'index': 36
				}
				,
				{
					'name': 'Mtb46',
					'index': 37
				}
				,
				{
					'name': 'Mtb47',
					'index': 38
				}
				,
				{
					'name': 'Mtb49',
					'index': 39
				}
				,
				{
					'name': 'Mtb50',
					'index': 40
				}
				,
				{
					'name': 'Mtb51',
					'index': 41
				}
				,
				{
					'name': 'Mtb52',
					'index': 42
				}
				,
				{
					'name': 'Mtb53',
					'index': 43
				}
				,
				{
					'name': 'Mtb54',
					'index': 44
				}
				,
				{
					'name': 'Mtb55',
					'index': 45
				}
				,
				{
					'name': 'Mtb57',
					'index': 46
				}
				,
				{
					'name': 'Mtb59',
					'index': 47
				}
				,
				{
					'name': 'Mtb68',
					'index': 48
				}
				,
				{
					'name': 'Mtb69',
					'index': 49
				}
				,
				{
					'name': 'Mtb70',
					'index': 50
				}
				,
				{
					'name': 'Mtb71',
					'index': 51
				}
				,
				{
					'name': 'Mtb72',
					'index': 52
				}
				,
				{
					'name': 'Mtb74',
					'index': 53
				}
				,
				{
					'name': 'Mtb75',
					'index': 54
				}
				,
				{
					'name': 'Mtb77',
					'index': 55
				}
				,
				{
					'name': 'Mtb78',
					'index': 56
				}
				,
				{
					'name': 'Mtb80',
					'index': 57
				}
				,
				{
					'name': 'Mtb81',
					'index': 58
				}
				,
				{
					'name': 'Mtb82',
					'index': 59
				}
				,
				{
					'name': 'Mtb83',
					'index': 60
				}
				,
				{
					'name': 'Mtb84',
					'index': 61
				},
				{
					'name': 'Mtb85',
					'index': 62
				},
				{
					'name': 'Mtb86',
					'index': 63
				},
				{
					'name': 'Mtb87',
					'index': 64
				},
				{
					'name': 'Mtb88',
					'index': 65
				},
				{
					'name': 'Mtb90',
					'index': 66
				}

			]
		}
	},
	'associations': {
		'file' : [
		{	
			'condition' : 'Stimulated',
			'file' : 'asso.stim.MTB.csv'
		},
		{
			'condition' : 'Not Stimulated',
			'file' : 'asso.not_stim.MTB.csv'
		}

		]
	}	
} 
# config_Zeller = {
	
# 	project_name : 'Zeller',
# 	'association' : {


# 	}


# }

# A Snp document contains the snp_id, the location (chromosome and position), both possible alleles. It also contains a field genotype which is a dictionnary.
# For each project, and among each project, each condition, all genotypes are listed for the set of individuals with their labels. 
class Genotype_values(Document):
	condition_name = StringField()
	individuals = ListField(DictField(StringField()))

class Genotype(Document):
	project_name = StringField()
	geno0 = StringField()
	geno1 = StringField()
	geno2 = StringField()
	conditions = ListField(DocumentField(Genotype_values), db_field='conditions')

class Snp(Document):
	config_collection_name ='snps_genotypes'
	# _id = StringField()
	snp_id = StringField()
	snp_chr = StringField()
	snp_position = IntField()
	# mongo_id = snp_id
	genotype = ListField(DocumentField(Genotype), db_field = 'genotype')
	#Create indexes
	index_snp_id = Index().ascending('snp_id')
	index_snp_chr = Index().ascending('snp_chr').ascending('snp_position')



class Expression_condition(Document):
	condition = StringField()
	individuals = ListField(DictField(StringField()))

class Expression(Document):
	project_name = StringField()
	conditions = ListField(DocumentField(Expression_condition),db_field='conditions')

class Exon(Document):
	ensembl_exon_id = StringField()
	exon_start = IntField()
	exon_end = IntField()


# A gene contains the ensembl_gene_id and the external_id as well. It also contains the location of that gene (chromosome, position1 and position2) all exons 
# of a given gene are listed as a dictionnary where an exon is defined by its exon_id and its location. The field expression contains for each condition in each project
# The set of individuals are listed with their expression values for a given gene.  
class Gene(Document):
	config_collection_name = "genes_expression"
	external_gene_id = StringField()
	ensembl_gene_id = StringField()
	_id = ensembl_gene_id
	gene_chr = StringField()
	gene_start = IntField()
	gene_end = IntField()
	exons = ListField(DocumentField(Exon), db_field='exons')
	expression = ListField(DocumentField(Expression), db_field='expression')
	#Create indexes 
	index_ensembl_gene_id = Index().ascending('ensembl_gene_id')
	index_external_gene_id = Index().ascending('external_gene_id')
	index_chr_pos1_pos2 = Index().ascending('gene_chr').ascending('gene_start').ascending('gene_end')
	index_gene_chr = Index().ascending('gene_chr')
	index_position1 = Index().ascending('gene_start')
	index_position2 = Index().ascending('gene_end')




# An association is defined by the snp_id and the ensembl_gene_id + external_gene_id. The field association contains, for each condition in each projet
# the slope of the reqtl and the pvalue of the association in the form of a dictionnary
# class Association_condition(Document):
# 	condition = StringField()
# 	estimate = StringField()
# 	pval = StringField()
# 	type = StringField()

class Association_value(Document):
	project_name = StringField()
	condition = StringField()
	estimate = StringField()
	pval = StringField()
	type = StringField()

class Association(Document):
	config_collection_name = "associations"
	snp_id = StringField()
	ensembl_gene_id = StringField()
	external_gene_id = StringField()
	association = ListField(DocumentField(Association_value),db_field='association')
	#Create indexes
	index_snp_id = Index().ascending('snp_id')
	index_ensembl_gene_id = Index().ascending('ensembl_gene_id')
	index_snp_n_ensembl = Index().ascending('snp_id').ascending('ensembl_gene_id')


def add_associations(config):
	session = Session.connect('immunio')				
	# Feedback during execution 
	ass = []	
	print 'Starting to add Associations'
	first = True
	asso = Association()	
	for it in config['associations']['file']:
		print 'Starting to add Associations from file %s' % (it['file'])
		f = open(it['file'],'r')
		con = it['condition']
		# first = True
		for lines in f:
			line = lines.strip('\n').split(',')
			this_asso = {
				'snp_id':line[3],
				'ensembl_gene_id':line[0],
				'external_gene_id':line[1],
				'association':[],
				'estimate':line[6],
				'pval':line[7],
				'asso_type':line[8]
			}
			asso_con = Association_value()
			asso_con.condition = con
			asso_con.project_name = config['project_name']
			asso_con.estimate = this_asso['estimate']
			asso_con.pval = this_asso['pval']
			asso_con.type = this_asso['asso_type']
			a = Association(snp_id=this_asso['snp_id'],ensembl_gene_id=this_asso['ensembl_gene_id'],external_gene_id=this_asso['external_gene_id'], association=[asso_con])
			# asso_project.conditions.append(asso_con)
			# if first:
				# a_p = Association_value(project_name=config['project_name'],conditions = [])
				# a_p.conditions.append(asso_con)
				# a = Association(snp_id = this_asso['snp_id'], ensembl_gene_id=this_asso['ensembl_gene_id'], external_gene_id = this_asso['external_gene_id'], association=[a_p])
				# ass.append(a)
			# else:
				# found = False
				# for it in ass:
					# if it.ensembl_gene_id == this_asso['ensembl_gene_id'] and it.snp_id==this_asso['snp_id']:
						# it['association'].append(asso_con) 
						# found = True 
						# continue
				# if not found:
					# a_p = Association_value(project_name=config['project_name'],conditions = [])
					# a_p.conditions.append(asso_con)
					# a = Association(snp_id = this_asso['snp_id'], ensembl_gene_id=this_asso['ensembl_gene_id'], external_gene_id = this_asso['external_gene_id'], association=[a_p])
					# ass.append(a)
		# first = False
			associations = session.query(Association).filter({'snp_id':it['snp_id'], 'ensembl_gene_id':it['ensembl_gene_id']})
			if len(associations.all())==0:
				session.save(a)
			else:
				associations.append('association',asso_con).execute()
	print 'Done'
# def add_associations(config):
# 	session = Session.connect('immunio')				
# 	# Feedback during execution 
# 	ass = []	
# 	print 'Starting to add Associations'
# 	asso = Association()	
# 	asso_project = Association_value(project_name=config['project_name'],conditions = [])
# 	for it in config['associations']['file']:
# 		print 'Starting to add Associations from file %s' % (it['file'])
# 		f = open(it['file'],'r')
# 		con = it['condition']
# 		first = True
# 		for lines in f:
# 			line = lines.strip('\n').split(',')
# 			asso.snp_id=line[3]
# 			asso.ensembl_gene_id=line[0]
# 			asso.external_gene_id=line[1]
# 			asso.association=[]
# 			estimate = line[6]
# 			pval = line[7]
# 			asso_type = line[8]
# 			asso_con = Association_condition()
# 			asso_con.condition = con
# 			asso_con.estimate = estimate
# 			asso_con.pval = pval
# 			asso_con.type = asso_type
# 			asso_project.conditions.append(asso_con)
# 			associations = session.query(Association).filter({'snp_id':line[3], 'ensembl_gene_id':line[0]})		
# 			if not first:
# 				if len(associations.all())==0:
# 					asso.association.append(asso_project)
# 					session.save(asso)
# 				else:
# 					associations.append('association',asso).execute()
# 			else:
# 				first = False
# 	#Feedback during execution 
# 	print 'Done'



# This method allows to fill the mongodb with referenced snps. The field genotype is left empty as this step is just to build the reference. The file read 
# is a csv where the first line is the headers and the order of columns is snp_chr, position, snp_id, allele1 and allele2
def fill_snps(config):
	session = Session.connect('immunio')			
	f= open(config['reference']['snp_file'],"r")
	headers = False
	print 'Starting to add Snps from Reference to MongoDB browser'
	for lines in f:
		if not headers:
			line = lines.strip("\n").split("\t")
			chrom = line[0] 
			pos = line[1]
			snp_id = line[2]
			session.save(Snp(snp_id=snp_id, snp_chr=chrom, genotype = [], snp_position=int(pos)))
		headers = False
	print 'All Snps Added'

# This method fills the snps_genotypes collection with the genotypes of the individuals within a given project. The file should be specified in a config 
# structure 
def add_genotypes(config):
	session = Session.connect('immunio')				
	f = open(config['genotypes']['file'], "r")
	genotype = Genotype()
	header = True
	print 'Starting to Add genotypes for each Snp'
	for lines in f:
		if not header:
			line = lines.strip("\n").split(",")
			genotype.project_name = config['project_name']
			genotype.geno0 = line[3]
			genotype.geno1 = line[4]
			genotype.geno2 = line[5]
			genotype.conditions = []
			# genotype = dict([["project_name",config['project_name']],["conditions",[]]])
			individuals = [] 
			for it in config['genotypes']['fields']['individuals']:
				individuals.append(dict([["individual",it['name']],["genotype_code",line[it['index']]]]))	
			count = 0	
			for con in config['genotypes']['fields']['conditions']:
				genotype_condition = Genotype_values()
				genotype_condition.condition_name = con
				genotype_condition.individuals = individuals
				genotype.conditions.append(genotype_condition)
				# genotype['conditions'].append(dict([["condition",condition[count]],["values",individuals]])) 
				count += 1
			query = session.query(Snp).filter(Snp.snp_id == line[config['genotypes']['fields']['snp_id']])
			# Verifier la correspondance entre les genotypes et les alleles de references
			if len(query.raw_output().all())==0:
				snp_id = line[0]
				chrom = line[1]
				pos = line[2]
				session.save(Snp(snp_id=snp_id, snp_chr=chrom, genotype = [genotype], snp_position=int(pos)))
			else:
				query.append('genotype', genotype).execute()
		else:
			header = False 
	print 'Done'


# This method fills the expressions collection with the expression for each gene referenced in the latter collection. For each genes, 
# the method loops on all individuals and a dictionnary containing the condition for the given project in the form of a dictionnary. 

# This method builds the reference for all genes. A coordinate file needs to be specified into the config_ref and all genes in this file 
# Are loaded in the expressions collection. To add all exons for a given gene, we need to query the sqlite3 database which contains all exons 
# linked with a corresponding ensembl_gene_id
def fill_genes(config):
	session = Session.connect('immunio')			
	f = open(config['reference']['gene_file'],"r")
	headers = True
	conn = sqlite3.connect(DB)
	conn.row_factory = sqlite3.Row
	db = conn.cursor()	
	# Feedback during execution 
	print 'Starting to add Genes from Reference to MongoDB browser'
	for lines in f:
		if not headers:
			line = lines.strip("\n").split(",")
			ensembl_gene_id = line[0] 
			external_gene_id = line[1]
			gene_chr = line[2]
			start = line[3]
			end = line[4]
			qryExons = 'SELECT * FROM Exons WHERE ensembl_gene_id==\'%s\''%(ensembl_gene_id)
			rowExon = db.execute(qryExons)
			exons = [dict(i) for i in rowExon]
			exons_out = []	
 			ge = Gene(_id = ensembl_gene_id, ensembl_gene_id=ensembl_gene_id, external_gene_id=external_gene_id, gene_chr=gene_chr, gene_start=int(start), gene_end=int(end), exons=[], expression= [])		
			for ex in exons:
				exon = Exon()
				exon.ensembl_exon_id = ex['ensembl_exon_id']
				exon.exon_start = int(ex['exon_start'])
				exon.exon_end = int(ex['exon_end'])
 				ge.exons.append(exon)
			session.save(ge)
		headers = False	
	conn.commit()
	conn.close()	
	# Feedback during execution 
	print 'All Genes Added'

def add_expressions(config):
	session = Session.connect('immunio')
	# Feedback during execution
	exprs = []
	print 'Starting to Add expressions for each individuals per Gene Tested'
	for f in config['expressions']['file']:
		expressions = open(f['name'],'r')
		header = True
		for ex in expressions:
			if not header:				
				expression = Expression(project_name=config['project_name'],conditions=[])
				individuals=[]
				line = ex.strip('\n').split(',')	
				ensembl = line[config['expressions']['fields']['ensembl_gene_id']]
				for it in config['expressions']['fields']['individuals']:
					individuals.append(dict([["individual",it['name']],["value",line[it['index']]]]))
				expr_con = Expression_condition(condition = f['condition'],individuals = individuals)
				found = False
				for item in exprs:
					if item['ensembl']==ensembl:
						item['expression'].conditions.append(expr_con)
						found = True
						break
				if not found:
					expression.conditions.append(expr_con)
					exprs.append(
						{
							'ensembl':ensembl,
							'expression':expression
						})
			else:
				header = False
		expressions.close()	
	for items in exprs:
		query = session.query(Gene).filter(Gene.ensembl_gene_id == items['ensembl'])
		query.append('expression', items['expression']).execute()	
		# query = session.query(Gene).filter(Gene.ensembl_gene_id == line[config['expressions']['fields']['ensembl_gene_id']])		
		# query.append('expression', expression).execute()		
	# Feedback during execution 
	print 'Done'

# This method goes through a dumped file from the sqlite3 db and adds every couple of snp/gene to the collection. If such association is absent from the collection 
# the association is added and the field association is appended with the specific values for this project/condition else, the values are simply appended to the 
# existing association with information about the specific project. 
if __name__=='__main__':
	# fill_snps(config_ref)
	# fill_genes(config_ref)
	# add_genotypes(config_MTB)
	# add_expressions(config_MTB)
	add_associations(config_MTB)



