from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import argparse, os, re
# 建立中间产物文件夹
# 输入blast_db文件，构建本地blast_db数据库
# 调用本地blast,采用input的方式个性化选择blast参数，解析blast输出文件
# 调用软件mafft，把query构建成，hmm文件.
# 修改基因组序列id,从output文件中提取序列id，调用Hmmer，解析Hmmer结果，
# 取两者交集，获得序列
class get_message():
    def __init__(self):
        par = argparse.ArgumentParser()
        par.add_argument('-query', '--query')  # 目标query
        par.add_argument('-f', '--file')  # 基因组文件
        arg = par.parse_args()
        self.query = arg.query
        self.file = arg.file
        os.mkdir('result')
        os.mkdir('result/midfile')
    def bio(self,path,way):
        if way == 'prot':
            t=[]
            for i in SeqIO.parse(path,'fasta'):
                a=Seq(str(i.seq)).upper()
                a=SeqRecord(a,id=i.id,description='')
                t.append(a)
            path='result/midfile/'+str(path).split('/')[-1]

            SeqIO.write(t,path,'fasta')
            return path
        else:
            t=[]
            for i in SeqIO.parse(path,'fasta'):
                try:
                    a=Seq(str(i.seq)).upper().translate()
                    a=SeqRecord(a,id=i.id,description='')
                    t.append(a)
                except:
                    print(i.id,'翻译出错')
                    break
            if len(t) != 0:
                path = 'result/midfile/' + str(path).split('/')[-1]
                SeqIO.write(t, path, 'fasta')
            return path
    def data_pre(self):
        file_type=input('序列类型(gene/prot):')
        if file_type=='prot':
            self.query=self.bio(self.query,'prot')
            self.file=self.bio(self.file,'prot')
        else:
            self.query=self.bio(self.query,'cds')
            self.file=self.bio(self.file,'cds')
    def blast_pipeline(self):
        print('---------------------------running blast----------------------------------')
        e = input('设定e值：')
        ident = int(input('最小序列一致性:'))
        try:
            os.system(f"makeblastdb -in {self.query} -out result/midfile/db -dbtype prot")
        except:
            print('构建本地数据库出错')
        # 调用本地blast
        try:
            outfile = 'result/midfile/blast_out.file'
            os.system(
                f'blastp -query {self.file} -db result/midfile/db -out {outfile} -outfmt "6 qseqid qstart qend sacc score pident " -evalue {e}')
        except:
            print('blast出错')
        df = pd.read_table(outfile, names=['条带编号', '起始位点', '终止位点', '酶编号', '得分', '一致性'])
        self.blast_name = list(set([i for i in df[df['一致性'] > ident]['条带编号']]))
        print('blast 条数：', len(self.blast_name))
        print('                                                 ')
    def Hmmer_pipeline(self):
        # 构建本地hmm文件，改序列名字，
        print('---------------------------running Hmmsearch----------------------------------')
        os.system(f'linsi --anysymbol {self.query} > result/midfile/mafft.fas')
        #a = []
        #for i in SeqIO.parse('result/midfile/mafft.fas', 'fasta'):
            #p = SeqRecord(i.seq, id=i.id, description='')
            #a.append(p)
        #SeqIO.write(a, 'result/midfile/hmm.sto', 'stockholm')
        SeqIO.convert('result/midfile/mafft.fas', 'fasta','result/midfile/hmm.sto', 'stockholm')


        
        # 构建本地hmm索引
        os.system(f'hmmbuild result/midfile/query.hmm result/midfile/hmm.sto')
        # 改序列名字
        a = []
        for i in SeqIO.parse(self.file, 'fasta'):
            p = SeqRecord(i.seq, id='YR_' + i.id, description='')
            a.append(p)
        SeqIO.write(a, 'result/midfile/hmm.fas', 'fasta')

        os.system(f'hmmsearch result/midfile/query.hmm result/midfile/hmm.fas > result/midfile/hmmsearch.txt')

        # 提取序列信息
        with open('result/midfile/hmmsearch.txt', 'r') as f:
            a = re.findall(r'YR_.*? ', f.read())
            self.hmm_name = list(set([i.strip() for i in a]))
            self.hmm_name = [i.replace('YR_', '') for i in self.hmm_name]

        print('hmmseach 条数：', len(self.hmm_name))
    def main(self):
        self.data_pre()
        self.blast_pipeline()
        self.Hmmer_pipeline()
       
        final_seq = self.sec_selcted_seq = [i for i in self.hmm_name if i in self.blast_name]

        a = {}
        for i in SeqIO.parse(self.file, 'fasta'):
            if str(i.id) in final_seq:
                a[str(i.id)]=str(i.seq)
        with open('result/final.fas','w') as f:
            for name in a:
                f.write('>'+name+'\n')
                f.write(a[name]+'\n')
        print('序列已保存为result/final.fas')
yr = get_message()
yr.main()
