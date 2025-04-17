import matplotlib.pyplot as plt
import math
import numpy as np

def emulate(states, reaction_list, render_list=None, accuracy=1000, finish=5e-5, max_time=2000, max_render_points=100):
	if render_list == None:
		render_list = list(states.keys())
	time_states = [0]
	t = 0
	dt = 1/accuracy
	print('EMULATION STARTED')
	while t<max_time:
		if int(t) % 10 == 0 and abs(t-int(t)) < dt:
			print("time:", t)
		sub_deltas = {}
		for key in states.keys():
			sub_deltas[key]=0.0

		for i in range(len(reaction_list)):
			delta = reaction_list[i][2] * dt
			for name in reaction_list[i][0]:
				delta *= states[name][-1]

			for name in reaction_list[i][0]:
				sub_deltas[name] -= delta
			for name in reaction_list[i][1]:
				sub_deltas[name] += delta

		for key in states.keys():
			states[key].append(states[key][-1] + sub_deltas[key])

		t += dt
		time_states.append(t)

		m = 0
		for value in sub_deltas.values():
			m = max(m, abs(value))
		if m < finish * dt:
			break

	print("EMULATION FINISHED\n total time:", t)

	scale = math.ceil(len(time_states)/max_render_points)

	for key in render_list:
		plt.plot(time_states[::scale], states[key][::scale], label=key)

	plt.xticks([])
	plt.yticks([0,1])

	plt.legend()
	plt.show()


if __name__ == '__main__':

	Mystates = {'RH':[1],
			  'Br2':[1],
			  'R1':[0],
			  'R2':[0],
			  'Br':[0],
			  'R1Br':[0],
			  'R2Br':[0],
			  'HBr':[0],
			  'R1R2':[0],
	}

	Myreaction_list = [[['Br2'], ['Br','Br'], 0.01],
					 [['Br', 'Br'], ['Br2'], 0.0],
					 [['Br', 'RH'], ['HBr', 'R1'], 5.0],
					 [['Br', 'RH'], ['HBr', 'R2'], 50.0],
					 [['R1', 'Br2'], ['R1Br', 'Br'], 0.3],
					 [['R2', 'Br2'], ['R2Br', 'Br'], 3.0],
					 [['HBr', 'R1'], ['Br', 'RH'], 1.0],
					 [['HBr', 'R2'], ['Br', 'RH'], 1.0],
					 [['R1Br', 'Br'], ['R1', 'Br2'], 3.0],
					 [['R2Br', 'Br'], ['R2', 'Br2'], 30.0],
					 [['R1', 'Br'], ['R1Br'], 100.0],
					 [['R2', 'Br'], ['R2Br'], 100.0],
					 # [['R1', 'R2'], ['R1R2'], 50.0],
	]

	emulate(Mystates, Myreaction_list, render_list=['R1Br', 'R2Br'], max_time=650)


