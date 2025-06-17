import matplotlib.pyplot as plt
from prettytable import PrettyTable
import numpy as np

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['font.size'] = 16
plt.rcParams['axes.unicode_minus'] = False

class Reation:
    """基元反应类"""

    def __init__(self,rectants:list, products:list, rate_constant,name='defalut'):
        """
        reactants: 一个包含所有反应物的列表,反应物用字符串表示
        products: 一个包含所有产物的列表,产物用字符串表示
        """
        self.reactants = rectants
        self.products = products
        self.rate_constant = rate_constant
        self.name = name

class DynamicModel:
    """总动力学模型"""

    def __init__(self,reation_list:list[Reation]):
        """reationlist:该动力学模型所包含的所有基元反应组成的列表"""
        if type(reation_list)!= list or len(reation_list)==0 or type(reation_list[0])!=Reation:
            raise AttributeError('invalid input')
        self.reation_list = reation_list
        self.pool = dict()
        self.initialization()
        self.logs = dict()
        for key in self.pool:
            self.logs[key] = [self.pool[key]]
        self.frozen_list = []
    
    def initialization(self,ini_con=dict()):
        """对所有底物浓度初始化,ini_con为一个储存初始浓度的字典.
        未在ini_con中给出浓度的反应物浓度将被初始化为1.0
        未在ini_con中给出浓度的产物浓度将被初始化为0.0
        """
        for reation in self.reation_list:
            for item in reation.reactants:
                self.pool[item]=1.0
            for item in reation.products:
                self.pool[item]=0.0
        for key in ini_con.keys():
            if key not in self.pool.keys():
                raise AttributeError(f"invalid input: {key}")
            self.pool[key] = ini_con[key]
            self.logs[key] = [ini_con[key]]

    def forzen(self, forzen_spices:list[str]):
        """输入一个包含物种名称字符串的列表,其中所有的物种的浓度在模拟过程中不会改变"""
        for key in forzen_spices:
            if key not in self.pool.keys():
                raise AttributeError(f"invalid input: {key}")
            self.frozen_list.append(key)

    def simulate(self, dt=0.01, stop_condition=1e-5, max_step=1e5):
        """
        启动动力学模拟的方法
        输入：dt : float ,模拟的时间单位
        stop_condition : int . 当反应速率最大的反应的速率小于该值时，认为模型收敛，模拟提前终止
        max_step : int , 当模拟步数达到该值时，停止模拟。
        输出：
        self.logs , 一个保存所有物种浓度-时间信息的字典
        """
        print('simulate start.')
        step = 1
        derta_max = stop_condition

        log_interval = max(1, int(max_step // 1e5)) if max_step > 1e5 else 1
        self.logs['time'] = [0.0]
        frozen_set = set(self.frozen_list)

        initial_length = len(self.logs[next(iter(self.pool))])
        need_initial_record = initial_length == 1
        while derta_max >= stop_condition and step <= max_step:
            derta_max = 0

            delta_pool = {s: 0.0 for s in self.pool}
            for reaction in self.reation_list:
                rate = reaction.rate_constant * dt
                for item in reaction.reactants:
                    rate *= self.pool[item]
                for item in reaction.reactants:
                    if item not in frozen_set:
                        delta_pool[item] -= rate
                for item in reaction.products:
                    if item not in frozen_set:
                        delta_pool[item] += rate
                derta_max = max(abs(rate), derta_max)

            for k, v in delta_pool.items():
                self.pool[k] += v

            if step % log_interval == 0 or step == max_step:
                for species in self.pool:
                    self.logs[species].append(self.pool[species])
                self.logs['time'].append(step * dt)
            step += 1

        if len(self.logs['time']) == 1:
            for species in self.pool:
                self.logs[species].append(self.pool[species])
            self.logs['time'].append(0.0)
        print(f'simulate finished. Cost {step-1} steps.')
        return self.logs
    
    def render(self, spices: list, scale_threshold: float = 1e2, name='浓度-时间图像', auto_log=True):
        """
        将模拟结果可视化
        输入： spices : list， 包含想要绘制图像的物种的名称的列表
        scale_threshold ：float , 当物种浓度的最大最小值之间的比值达到该值时，采用浓度的对数-时间图像
        auto_log : bool , 值为True时自动检测是否需要取对数
        name: str ， 图表名称
        """
        for species in spices:
            if species not in self.logs:
                raise AttributeError(f'Invalid species: {species}')

        time = self.logs['time']
        plt.figure(figsize=(12, 8))

        transform_type = 'linear'
        if auto_log and spices:
            conc_avgs = []
            for species in spices:
                conc = np.nan_to_num(self.logs[species], nan=0.0)
                non_zero_conc = conc[conc > 0]
                if len(non_zero_conc) > 0:
                    conc_avgs.append(np.mean(non_zero_conc))
                else:
                    conc_avgs.append(0.0)

            nonzero_avgs = [avg for avg in conc_avgs if avg > 0]
            if nonzero_avgs:
                max_conc = max(nonzero_avgs)
                min_conc = min(nonzero_avgs)
                if (max_conc / min_conc) >= scale_threshold:
                    transform_type = 'log'
            
        ax = plt.subplot(1, 1, 1)
        for species in spices:
            conc = np.nan_to_num(self.logs[species], nan=0.0)
            step = max(1, len(conc) // 1000)

            if transform_type == 'log':
                data = np.log10(np.where(conc > 0, conc, 1e-10))
                ylabel = 'log10(浓度)'
            else:
                data = conc
                ylabel = '浓度'

            ax.plot(time[::step], data[::step],
                linestyle='-',
                label=f'{species} ({ylabel})')

        ax.set_ylabel(ylabel)
        ax.set_xlabel('时间')
        ax.legend(loc='best')
        ax.grid(True)
        plt.title(name)
        plt.tight_layout()
        plt.show()

    def summary(self):
        table = PrettyTable()
        table.field_names = ["Reaction", "Rate Constant", "Initial Concentrations"]

        for reaction in self.reation_list:
            reactants_str = " + ".join(reaction.reactants)
            products_str = " + ".join(reaction.products)
            reaction_str = f"{reactants_str} → {products_str}"

            init_conc = []
            for r in reaction.reactants:
                init_conc.append(f"{r}: {self.logs[r][0]:.4f}")
            init_conc_str = ", ".join(init_conc)

            table.add_row([reaction_str, reaction.rate_constant, init_conc_str])

        print("Reaction System Summary:")
        print(table)

        print('These spices are forzen:')
        for i in self.frozen_list: 
            print(i,end = ',')
            print()
