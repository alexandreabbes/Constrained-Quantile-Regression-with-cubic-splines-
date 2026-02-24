#!/usr/bin/env python

__author__ = "Alexandre Abbes"
__copyright__ = "Copyright 2026, Alexandre Adel Abbes"
__credits__ = ["Alexandre Abbes"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Alexandre Abbes"
__email__ = "alexandre.abbes@proton.me"


import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import warnings
warnings.filterwarnings('ignore')

#Cette version inclut les regressions sous contraintes pour des splines de tous les degrés de 1 (affines) à 4 (quartiques).
#Les  contraintes sont choisies globales ou par régions pour des degrés de 1 à 3 maximum. Les contraintes de degré 3 sont abordées à titre expérimental pour les splines de degé 3 et 4, les contraintes de degré 4 ne sont pas abordées.

# Import direct de vos fonctions
try:
    from SplineCubicQuantBspkn3_V3_der3 import SplineCubicQuantBspkn3
    print("✓ Module spline cubique chargé")
except ImportError as e:
    print(f"✗ Erreur chargement spline cubique: {e}")
    SplineCubicQuantBspkn3 = None

try:
    from SplineQuarticQuantBspkn4_V3 import SplineQuantBspkn4_with_der3 as SplineQuantBspkn4
    print("✓ Module spline quartique chargé")
except ImportError as e:
    print(f"✗ Erreur chargement spline quartique: {e}")
    SplineQuantBspkn4 = None

try:
    from SplineQuadQuantBspkn2_V1 import SplineQuadraticQuantile
    print("✓ Module spline quadratique chargé")
except ImportError as e:
    print(f"✗ Erreur chargement spline quadratique: {e}")
    SplineQuadraticQuantile = None

try:
    from SplineLinQuantBspkn1_V1 import SplineLinearQuantile
    print("✓ Module spline linéaire chargé")
except ImportError as e:
    print(f"✗ Erreur chargement spline linéaire: {e}")
    SplineLinearQuantile = None

    
class QuantileSplineApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Régression Quantile avec Splines Contraintes")
        self.root.geometry("1400x800")
        
        # Données
        self.xtab = None
        self.ytab = None
        self.knots = None
        self.manual_knots_mode = False
        self.temp_knots = []
        self.knot_lines = []
        self.spline_result = None
        self.degree = 3
        self.constraints_mode = "uniform"
        self.regions = []
        self.region_patches = []
        
        # Solveurs disponibles
        self.solvers = [
            'GUROBI',
            'MOSEK', 
            'CLARABEL',
            'ECOS',
            'SCS',
            'CVXOPT'
        ]
        
        self.setup_ui()
        
    def setup_ui(self):
        """Configuration de l'interface"""
        
        # Frame principal
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configuration des poids
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=3)
        main_frame.rowconfigure(0, weight=1)
        
        # ============ NOTEBOOK avec onglets (à gauche) ============
        notebook = ttk.Notebook(main_frame)
        notebook.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)
        
        # Configuration du poids pour que le notebook s'étende
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(0, weight=1)
        
        # ============ Onglet 1: Données et paramètres spline ============
        tab1 = ttk.Frame(notebook)
        notebook.add(tab1, text="1. Données & Spline")
        
        # ============ Onglet 2: Contraintes ============
        tab2 = ttk.Frame(notebook)
        notebook.add(tab2, text="2. Contraintes")
        
        # ============ Onglet 3: Exécution ============
        tab3 = ttk.Frame(notebook)
        notebook.add(tab3, text="3. Exécution")
        
        # ------------------------------------------------------------
        # CONTENU DE L'ONGLET 1
        # ------------------------------------------------------------
        
        # Section Import des données
        import_frame = ttk.LabelFrame(tab1, text="Import des données", padding="10")
        import_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=5, padx=5)
        
        ttk.Button(import_frame, text="Charger CSV", command=self.load_csv).grid(row=0, column=0, padx=2)
        ttk.Button(import_frame, text="Charger Excel", command=self.load_excel).grid(row=0, column=1, padx=2)
        ttk.Button(import_frame, text="Données test", command=self.generate_test_data).grid(row=0, column=2, padx=2)
        
        self.data_info = tk.StringVar(value="Aucune donnée")
        ttk.Label(import_frame, textvariable=self.data_info).grid(row=1, column=0, columnspan=3, pady=5)
        
        # Section Paramètres spline
        spline_frame = ttk.LabelFrame(tab1, text="Paramètres spline", padding="10")
        spline_frame.grid(row=2, column=0, sticky=(tk.W, tk.E), pady=5, padx=5)
        
        # Degré

        ttk.Label(spline_frame, text="Degré:").grid(row=0, column=0, sticky=tk.W)
        self.degree_var = tk.IntVar(value=3)
        ttk.Radiobutton(spline_frame, text="1 (Linéaire)", variable=self.degree_var, value=1, command=self.update_degree).grid(row=0, column=1, sticky=tk.W)
        ttk.Radiobutton(spline_frame, text="2 (Quadratique)", variable=self.degree_var,value=2, command=self.update_degree).grid(row=0, column=2, sticky=tk.W)
        ttk.Radiobutton(spline_frame, text="3 (Cubique)", variable=self.degree_var, value=3, command=self.update_degree).grid(row=1, column=1, sticky=tk.W)
        ttk.Radiobutton(spline_frame, text="4 (Quartique)", variable=self.degree_var, value=4, command=self.update_degree).grid(row=1, column=2, sticky=tk.W)

       
        # Mode nœuds
        ttk.Label(spline_frame, text="Nœuds:").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.knot_mode = tk.StringVar(value="auto")
        ttk.Radiobutton(spline_frame, text="Auto (quantiles)", variable=self.knot_mode, value="auto", command=self.knot_mode_changed).grid(row=2, column=1, columnspan=2, sticky=tk.W)
        ttk.Radiobutton(spline_frame, text="Manuel (clics)", variable=self.knot_mode, value="manual", command=self.knot_mode_changed).grid(row=3, column=1, columnspan=2, sticky=tk.W)
        # Nombre de nœuds (mode auto)
        ttk.Label(spline_frame, text="Nombre:").grid(row=4, column=0, sticky=tk.W, pady=5)
        self.knot_count = tk.IntVar(value=11)
        self.knot_spinbox = ttk.Spinbox(spline_frame, from_=4, to=30, textvariable=self.knot_count, width=5)
        self.knot_spinbox.grid(row=4, column=1, sticky=tk.W)
        
        # Boutons nœuds
        knot_btn_frame = ttk.Frame(spline_frame)
        knot_btn_frame.grid(row=5, column=0, columnspan=3, pady=5)
        
        self.define_knots_btn = ttk.Button(knot_btn_frame, text="Définir nœuds", command=self.define_knots)
        self.define_knots_btn.pack(side=tk.LEFT, padx=2)
        
        self.validate_knots_btn = ttk.Button(knot_btn_frame, text="Valider nœuds", command=self.validate_manual_knots)
        self.validate_knots_btn.pack(side=tk.LEFT, padx=2)
        self.validate_knots_btn.config(state='disabled')
        
        ttk.Button(knot_btn_frame, text="Effacer nœuds", command=self.clear_knots).pack(side=tk.LEFT, padx=2)
        test_frame = ttk.Frame(import_frame)
        test_frame.grid(row=2, column=0, columnspan=3, pady=5, sticky=(tk.W, tk.E))

        ttk.Label(test_frame, text="Fonction test (x):").grid(row=0, column=0, sticky=tk.W)

                # Variable pour la fonction test
        self.test_function_var = tk.StringVar(value="2*x + 0.2*sin(10*pi*x) + 0.05*randn(n)")

                # Champ de saisie
        test_entry = ttk.Entry(test_frame, textvariable=self.test_function_var, width=40)
        test_entry.grid(row=1, column=0, padx=2, pady=2, sticky=(tk.W, tk.E))

                # Bouton pour générer
        ttk.Button(test_frame, text="Générer test",command=self.generate_custom_test).grid(row=1, column=1, padx=2)

                # Info nombre de points
        ttk.Label(test_frame, text="Nb points:").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.test_npoints_var = tk.IntVar(value=200)
        test_npoints_spin = ttk.Spinbox(test_frame, from_=10, to=1000,
                                        textvariable=self.test_npoints_var, width=8)
        test_npoints_spin.grid(row=2, column=1, sticky=tk.W, pady=2)

        # ------------------------------------------------------------
        # CONTENU DE L'ONGLET 2
        # ------------------------------------------------------------
        
        # Section Contraintes
        constraint_frame = ttk.LabelFrame(tab2, text="Contraintes de forme", padding="10")
        constraint_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N), pady=5, padx=5)
        
        # Mode contraintes
        self.constraint_mode_var = tk.StringVar(value="uniform")
        ttk.Radiobutton(constraint_frame, text="Uniformes", variable=self.constraint_mode_var, 
                       value="uniform", command=self.toggle_constraint_mode).grid(row=0, column=0, sticky=tk.W)
        ttk.Radiobutton(constraint_frame, text="Par région", variable=self.constraint_mode_var, 
                       value="region", command=self.toggle_constraint_mode).grid(row=0, column=1, sticky=tk.W)
        
        # Frame contraintes uniformes
        self.uniform_frame = ttk.Frame(constraint_frame)
        self.uniform_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        
        # Monotonie
        ttk.Label(self.uniform_frame, text="Monotonie:").grid(row=0, column=0, sticky=tk.W)
        self.monot_var = tk.IntVar(value=0)
        ttk.Radiobutton(self.uniform_frame, text="↗", variable=self.monot_var, value=1).grid(row=0, column=1)
        ttk.Radiobutton(self.uniform_frame, text="↘", variable=self.monot_var, value=-1).grid(row=0, column=2)
        ttk.Radiobutton(self.uniform_frame, text="✗", variable=self.monot_var, value=0).grid(row=0, column=3)
        
        # Convexité
        self.convex_frame = ttk.Frame(self.uniform_frame)

        self.convex_frame.grid(row=1, column=0, columnspan=4, sticky=(tk.W, tk.E), pady=5)
        ttk.Label(self.convex_frame, text="Convexité:").grid(row=0, column=0, sticky=tk.W)
        self.convex_var = tk.IntVar(value=0)
        ttk.Radiobutton(self.convex_frame, text="∪", variable=self.convex_var, value=1).grid(row=0, column=1)
        ttk.Radiobutton(self.convex_frame, text="∩", variable=self.convex_var, value=-1).grid(row=0, column=2)
        ttk.Radiobutton(self.convex_frame, text="✗", variable=self.convex_var, value=0).grid(row=0, column=3)        
        
        # Dérivée 3
        self.der3_frame = ttk.Frame(self.uniform_frame)
        self.der3_frame.grid(row=2, column=0, columnspan=4, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(self.der3_frame, text="Dérivée 3e:").grid(row=0, column=0, sticky=tk.W)
        self.der3_var = tk.IntVar(value=0)
        ttk.Radiobutton(self.der3_frame, text="+", variable=self.der3_var, value=1).grid(row=0, column=1)
        ttk.Radiobutton(self.der3_frame, text="-", variable=self.der3_var, value=-1).grid(row=0, column=2)
        ttk.Radiobutton(self.der3_frame, text="✗", variable=self.der3_var, value=0).grid(row=0, column=3)
        
        # Frame contraintes par région
        self.region_frame = ttk.Frame(constraint_frame)
        self.region_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        self.region_frame.grid_remove()
        
        ttk.Button(self.region_frame, text="Ajouter région", command=self.start_region_selection).pack(pady=2)
        ttk.Button(self.region_frame, text="Effacer régions", command=self.clear_regions).pack(pady=2)
        
        self.region_listbox = tk.Listbox(self.region_frame, height=4, width=35)
        self.region_listbox.pack(pady=5, fill=tk.X)
        
        # ------------------------------------------------------------
        # CONTENU DE L'ONGLET 3
        # ------------------------------------------------------------
        
        # Section Solveur
        solver_frame = ttk.LabelFrame(tab3, text="Solveur", padding="10")
        solver_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=5, padx=5)
        
        ttk.Label(solver_frame, text="Solveur CVXPY:").grid(row=0, column=0, sticky=tk.W)
        self.solver_var = tk.StringVar(value='GUROBI')
        solver_combo = ttk.Combobox(solver_frame, textvariable=self.solver_var, 
                                   values=self.solvers, state='readonly', width=15)
        solver_combo.grid(row=0, column=1, padx=5, pady=5)
        
        # Section Exécution
        exec_frame = ttk.LabelFrame(tab3, text="Exécution", padding="10")
        exec_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5, padx=5)
        
        # Tau
        ttk.Label(exec_frame, text="τ:").grid(row=0, column=0)
        self.tau_var = tk.DoubleVar(value=0.5)
        tau_scale = ttk.Scale(exec_frame, from_=0.01, to=0.99, variable=self.tau_var,
                            orient=tk.HORIZONTAL, length=150)
        tau_scale.grid(row=0, column=1, padx=5)
        self.tau_label = ttk.Label(exec_frame, text="0.50")
        self.tau_label.grid(row=0, column=2)
        tau_scale.configure(command=lambda x: self.tau_label.configure(text=f"{float(x):.2f}"))
        
        # Boutons
        btn_frame = ttk.Frame(exec_frame)
        btn_frame.grid(row=1, column=0, columnspan=3, pady=10)
        
        ttk.Button(btn_frame, text="Lancer", command=self.run_regression).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Exporter", command=self.export_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Effacer tout", command=self.clear_all).pack(side=tk.LEFT, padx=5)
        
        # Barre de progression
        self.progress = ttk.Progressbar(tab3, mode='indeterminate')
        self.progress.grid(row=2, column=0, sticky=(tk.W, tk.E), pady=10, padx=5)
        
        # Statut
        self.status_var = tk.StringVar(value="Prêt")
        ttk.Label(tab3, textvariable=self.status_var, foreground="blue").grid(row=3, column=0, pady=5, padx=5)
        

         # Cadre pour gestion des courbes
        curve_frame = ttk.LabelFrame(tab3, text="Gestion des courbes", padding="10")
        curve_frame.grid(row=2, column=0, sticky=(tk.W, tk.E), pady=5, padx=5)

         # Choix de la couleur
        color_frame = ttk.Frame(curve_frame)
        color_frame.grid(row=0, column=0, columnspan=2, pady=5, sticky=(tk.W, tk.E))

        ttk.Label(color_frame, text="Couleur:").grid(row=0, column=0, padx=2)
        self.color_var = tk.StringVar(value='blue')

        colors = ['blue', 'red', 'black', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
        color_combo = ttk.Combobox(color_frame, textvariable=self.color_var, 
                          values=colors, state='readonly', width=10)
        color_combo.grid(row=0, column=1, padx=2)

         # Bouton pour appliquer la couleur
        ttk.Button(color_frame, text="Appliquer", 
        command=self.apply_color).grid(row=0, column=2, padx=5)

         # Bouton pour effacer la dernière courbe
        curve_btn_frame = ttk.Frame(curve_frame)
        curve_btn_frame.grid(row=1, column=0, columnspan=2, pady=5)

        ttk.Button(curve_btn_frame, text="Effacer dernière courbe",  command=self.clear_last_curve).pack(side=tk.LEFT, padx=2)

        ttk.Button(curve_btn_frame, text="Effacer toutes les courbes",  command=self.clear_all_curves).pack(side=tk.LEFT, padx=2)

        # Compteur de courbes
        self.curve_count = 0
        self.curve_lines = []  # Pour stocker les lignes des courbes
        self.curve_count_var = tk.StringVar(value="Courbes: 0")
        ttk.Label(curve_frame, textvariable=self.curve_count_var).grid(row=2, column=0, columnspan=2, pady=5)
        # ============ Panneau de visualisation (droite) ============
        viz_frame = ttk.LabelFrame(main_frame, text="Visualisation", padding="10")
        viz_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)
        
        # Figure matplotlib
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_title('Régression quantile avec splines contraintes')
        self.ax.grid(True, alpha=0.3)
        self.ax.set_xlim([0, 1])
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=viz_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
        toolbar = NavigationToolbar2Tk(self.canvas, viz_frame)
        toolbar.update()
        
        # Événements souris
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        
        # Variables pour dessin
        self.press = None
        self.current_line = None
        
        # Initialisation
        self.toggle_constraint_mode()
        self.update_degree()
        self.knot_mode_changed()
        
        # Vérification des modules
        self.check_modules()
    
    def check_modules(self):
        """Vérifie que les modules sont disponibles"""
        if SplineCubicQuantBspkn3 is None:
            messagebox.showwarning("Attention", 
                "Module spline cubique non trouvé.\n"
                "Placez SplineCubicQuantBspkn3_V2_der3.py dans le même répertoire.")
        if SplineQuantBspkn4 is None:
            messagebox.showwarning("Attention", 
                "Module spline quartique avec derivée 3e non trouvé.\n"
                "Placez SplineQuarticQuantBspkn4_V3.py dans le même répertoire.")
       
            if SplineQuadraticQuantile is None:
                messagebox.showwarning("Attention", 
                                       "Module spline quadratique non trouvé.\n"
                                       "Placez SplineQuadQuantBspkn2_V1.py dans le même répertoire.")
            if SplineLinearQuantile is None:
                messagebox.showwarning("Attention", 
                                       "Module spline linéaire non trouvé.\n"
                                       "Placez SplineLinQuantBspkn1_V1.py dans le même répertoire.")

        # ============ Gestion des modes ============
    def toggle_constraint_mode(self):
        """Bascule mode contraintes"""
        mode = self.constraint_mode_var.get()
        if mode == "uniform":
            self.uniform_frame.grid()
            self.region_frame.grid_remove()
            self.constraints_mode = "uniform"
        else:
            self.uniform_frame.grid_remove()
            self.region_frame.grid()
            self.constraints_mode = "region"
    
    #def update_degree(self):
    #    """Met à jour l'interface selon le degré"""
    #    self.degree = self.degree_var.get()
    #    if self.degree == 3:
    #        self.der3_frame.grid()
    #    else:
    #        self.der3_frame.grid_remove()
    def update_degree(self):
        """Met à jour l'interface selon le degré"""
        self.degree = self.degree_var.get()

        # Afficher les contrôles appropriés selon le degré
        if self.degree >= 2:
            self.convex_frame.grid()
        else:
            self.convex_frame.grid_remove()
            self.convex_var.set(0)  # Remettre à zéro
    
            # Gérer la dérivée 3
        if self.degree >= 3:
                self.der3_frame.grid()
        else:
            self.der3_frame.grid_remove()
            self.der3_var.set(0)  # Remettre à zéro

        
                
    def knot_mode_changed(self):
        """Réagit au changement de mode de nœuds"""
        mode = self.knot_mode.get()
        if mode == "auto":
            self.knot_spinbox.config(state='normal')
            self.define_knots_btn.config(text="Définir nœuds", state='normal')
            self.validate_knots_btn.config(state='disabled')
            self.manual_knots_mode = False
        else:
            self.knot_spinbox.config(state='disabled')
            self.define_knots_btn.config(text="Démarrer saisie", state='normal')
            self.validate_knots_btn.config(state='normal')
            self.manual_knots_mode = False
    
    # ============ Import données ============
    def load_csv(self):
        filename = filedialog.askopenfilename(filetypes=[("CSV", "*.csv")])
        if filename:
            try:
                df = pd.read_csv(filename)
                self.process_dataframe(df)
            except Exception as e:
                messagebox.showerror("Erreur", str(e))
    
    def load_excel(self):
        filename = filedialog.askopenfilename(filetypes=[("Excel", "*.xlsx *.xls")])
        if filename:
            try:
                df = pd.read_excel(filename)
                self.process_dataframe(df)
            except Exception as e:
                messagebox.showerror("Erreur", str(e))
    
    def process_dataframe(self, df):
        """Traite le DataFrame"""
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        
        if len(numeric_cols) < 2:
            messagebox.showerror("Erreur", "Besoin de 2 colonnes numériques")
            return
        
        # Dialogue sélection
        dialog = tk.Toplevel(self.root)
        dialog.title("Sélection colonnes")
        dialog.geometry("300x200")
        dialog.transient(self.root)
        dialog.grab_set()
        
        ttk.Label(dialog, text="Colonne X:").grid(row=0, column=0, padx=5, pady=5)
        x_var = tk.StringVar(value=numeric_cols[0])
        ttk.Combobox(dialog, textvariable=x_var, values=numeric_cols, 
                    state="readonly").grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(dialog, text="Colonne Y:").grid(row=1, column=0, padx=5, pady=5)
        y_var = tk.StringVar(value=numeric_cols[1] if len(numeric_cols) > 1 else numeric_cols[0])
        ttk.Combobox(dialog, textvariable=y_var, values=numeric_cols,
                    state="readonly").grid(row=1, column=1, padx=5, pady=5)
        
    def confirm():
        try:
                self.xtab = df[x_var.get()].values.astype(float)
                self.ytab = df[y_var.get()].values.astype(float)
                
                # Normalisation [0,1]
                x_min, x_max = self.xtab.min(), self.xtab.max()
                if x_min != x_max:
                    self.xtab = (self.xtab - x_min) / (x_max - x_min)
                
                dialog.destroy()
                self.clear_all()
                self.plot_data()
                self.data_info.set(f"{len(self.xtab)} points")
                self.status_var.set("Données chargées")
        except Exception as e:
                 messagebox.showerror("Erreur", str(e))
        
        ttk.Button(dialog, text="Confirmer", command=confirm).grid(row=2, column=0, columnspan=2, pady=20)
    
    def generate_test_data(self):
                    """Données de test"""
                    np.random.seed(42)
                    n = 200
                    self.xtab = np.linspace(0, 1, n)
                    self.ytab = 2*self.xtab + 0.2*np.sin(10*np.pi*self.xtab) + 0.05*np.random.randn(n)
                    self.clear_all()
                    self.plot_data()
                    self.data_info.set(f"{n} points (test)")
                    self.status_var.set("Données test générées")
    
    # ============ Gestion nœuds ============
    def define_knots(self):
        """Définit les nœuds selon le mode"""
        if self.xtab is None:
            messagebox.showwarning("Attention", "Chargez d'abord des données")
            return
        
        if self.knot_mode.get() == "auto":
            # Mode automatique
            kn = self.knot_count.get()
            self.knots = np.quantile(self.xtab, np.linspace(0, 1, kn + 1))
            self.plot_knots()
            self.status_var.set(f"{len(self.knots)} nœuds auto (quantiles)")
        else:
            # Mode manuel: activer la saisie
            self.manual_knots_mode = True
            self.temp_knots = []
            self.status_var.set("Mode manuel: cliquez sur le graphique pour ajouter des nœuds")
            self.define_knots_btn.config(state='disabled')
    
    def validate_manual_knots(self):
        """Valide les nœuds saisis manuellement"""
        if len(self.temp_knots) < 2:
            messagebox.showwarning("Attention", "Ajoutez au moins 2 nœuds")
            return
        
        self.knots = np.array(sorted(self.temp_knots))
        self.manual_knots_mode = False
        self.define_knots_btn.config(state='normal')
        self.plot_knots()
        self.status_var.set(f"{len(self.knots)} nœuds manuels validés")
    
    def clear_knots(self):
        """Efface tous les nœuds"""
        self.knots = None
        self.temp_knots = []
        self.manual_knots_mode = False
        self.define_knots_btn.config(state='normal')
        
        # Supprimer les lignes
        for line in self.knot_lines:
            try:
                line.remove()
            except:
                pass
        self.knot_lines = []
        
        self.canvas.draw()
        self.status_var.set("Nœuds effacés")
    
    def plot_knots(self):
        """Affiche les nœuds"""
        # Effacer les anciens
        for line in self.knot_lines:
            try:
                line.remove()
            except:
                pass
        self.knot_lines = []
        
        if self.knots is not None:
            ylim = self.ax.get_ylim()
            y_pos = ylim[1] * 0.95
            
            for k in self.knots:
                line = self.ax.axvline(x=k, color='red', linestyle='--', alpha=0.5)
                self.knot_lines.append(line)
            
            # Marqueurs
            markers = self.ax.plot(self.knots, np.ones_like(self.knots) * y_pos, 
                                 'r|', markersize=10, label='Nœuds')
            self.knot_lines.extend(markers)
            
            self.ax.legend()
            self.canvas.draw()
    
    # ============ Gestion régions ============
    def start_region_selection(self):
        """Active sélection région"""
        self.constraint_mode_var.set("region")
        self.toggle_constraint_mode()
        self.status_var.set("Mode région: cliquez-glissez pour définir une zone")
    
    def clear_regions(self):
        """Efface toutes les régions"""
        self.regions = []
        self.region_listbox.delete(0, tk.END)
        
        for patch in self.region_patches:
            try:
                patch.remove()
            except:
                pass
        self.region_patches = []
        
        self.canvas.draw()
        self.status_var.set("Régions effacées")
    
    # ============ Gestion résultats ============
    def clear_results(self):
        """Efface les résultats"""
        self.spline_result = None
        self.plot_data()
        if self.knots is not None:
            self.plot_knots()
    
    def clear_all(self):
        """Efface tout"""
        self.clear_knots()
        self.clear_regions()
        self.clear_results()
        self.plot_data()
        self.status_var.set("Tout effacé")
    
    # ============ Visualisation ============
    def plot_data(self):
        """Affiche les données"""
        self.ax.clear()
        if self.xtab is not None and self.ytab is not None:
            self.ax.scatter(self.xtab, self.ytab, alpha=0.5, s=20, label='Données')
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_title('Régression quantile avec splines contraintes')
        self.ax.grid(True, alpha=0.3)
        self.ax.set_xlim([0, 1])
        self.ax.legend()
        self.canvas.draw()
    def apply_color(self):
        """Change la couleur pour les prochaines courbes"""
        # Juste un feedback, la couleur sera utilisée au prochain run
        self.status_var.set(f"Couleur changée: {self.color_var.get()}")

    def clear_last_curve(self):
     """Efface la dernière courbe tracée"""
     if self.curve_lines:
        last_line = self.curve_lines.pop()
        try:
            last_line.remove()
        except:
            pass
        self.curve_count = len(self.curve_lines)
        self.curve_count_var.set(f"Courbes: {self.curve_count}")
        self.canvas.draw()
        self.status_var.set("Dernière courbe effacée")

    def clear_all_curves(self):
        """Efface toutes les courbes sauf les données et les nœuds"""
        # Sauvegarder les éléments à conserver
        data_points = None
        knot_lines_copy = []
    
        # Trouver et supprimer uniquement les courbes de résultat
        for line in self.curve_lines:
            try:
                line.remove()
            except:
                pass
    
        self.curve_lines = []
        self.curve_count = 0
        self.curve_count_var.set("Courbes: 0")
    
        # Réafficher les données et les nœuds
        self.plot_data()
        if self.knots is not None:
            self.plot_knots()
    
        self.status_var.set("Toutes les courbes effacées")
        
    def plot_result(self):

        """Affiche le résultat"""
        if self.spline_result is not None:
            x_eval = np.linspace(0, 1, 500)
            y_eval = self.spline_result(x_eval)
        
            # Utiliser la couleur choisie
            color = self.color_var.get()
            line = self.ax.plot(x_eval, y_eval, color=color, linewidth=2,label=f'τ={self.tau_var.get():.2f}')
        
        # Stocker la ligne
        self.curve_lines.append(line[0])
        self.curve_count = len(self.curve_lines)
        self.curve_count_var.set(f"Courbes: {self.curve_count}")
        
        self.ax.legend()
        self.canvas.draw()
    
    # ============ Événements souris ============
    def on_click(self, event):
        if event.inaxes != self.ax or self.xtab is None:
            return
        
        # Mode ajout manuel de nœuds
        if self.manual_knots_mode:
            x = event.xdata
            if 0 <= x <= 1:
                self.temp_knots.append(x)
                
                # Afficher le nœud temporaire
                line = self.ax.axvline(x=x, color='orange', linestyle=':', alpha=0.7)
                self.knot_lines.append(line)
                self.canvas.draw()
                
                self.status_var.set(f"Nœud temporaire: {x:.3f} ({len(self.temp_knots)})")
        
        # Mode sélection région
        elif self.constraints_mode == "region" and self.knots is not None:
            self.press = event.xdata
            if self.current_line:
                self.current_line.remove()
            self.current_line = self.ax.axvline(x=self.press, color='green', alpha=0.5)
            self.canvas.draw()
    
    def on_motion(self, event):
        if self.press is None or event.inaxes != self.ax:
            return
        if self.constraints_mode == "region" and self.current_line:
            self.current_line.set_xdata([self.press, event.xdata])
            self.canvas.draw()
    
    def on_release(self, event):
        if self.press is None or event.inaxes != self.ax:
            return
        
        if self.constraints_mode == "region" and self.knots is not None:
            x1, x2 = self.press, event.xdata
            if x1 < x2 and 0 <= x1 <= 1 and 0 <= x2 <= 1:
                self.show_region_dialog(x1, x2)
            
            if self.current_line:
                self.current_line.remove()
                self.current_line = None
            self.press = None
            self.canvas.draw()

    def show_region_dialog(self, xmin, xmax):
        """Dialogue pour ajouter une région"""
        dialog = tk.Toplevel(self.root)
        dialog.title(f"Région [{xmin:.2f}, {xmax:.2f}]")
        dialog.geometry("350x300")
        dialog.transient(self.root)
        dialog.grab_set()
    
        row = 0
    
        # Monotonie (tous degrés)
        ttk.Label(dialog, text="Monotonie:").grid(row=row, column=0, padx=5, pady=5, sticky=tk.W)
        monot_var = tk.IntVar(value=0)
        ttk.Radiobutton(dialog, text="↗", variable=monot_var, value=1).grid(row=row, column=1)
        ttk.Radiobutton(dialog, text="↘", variable=monot_var, value=-1).grid(row=row, column=2)
        ttk.Radiobutton(dialog, text="✗", variable=monot_var, value=0).grid(row=row, column=3)
        row += 1
    
        # Convexité (degré >= 2)
        convex_var = tk.IntVar(value=0)
        # Der3 (deg>=3)
        der3_var = tk.IntVar(value=0)
        if self.degree >= 2:
            ttk.Label(dialog, text="Convexité:").grid(row=row, column=0, padx=5, pady=5, sticky=tk.W)
            ttk.Radiobutton(dialog, text="∪", variable=convex_var, value=1).grid(row=row, column=1)
            ttk.Radiobutton(dialog, text="∩", variable=convex_var, value=-1).grid(row=row, column=2)
            ttk.Radiobutton(dialog, text="✗", variable=convex_var, value=0).grid(row=row, column=3)
            row += 1
    
            # Dérivée 3 (degré >= 3)
            
        if self.degree >= 3:
                ttk.Label(dialog, text="Dérivée 3e:").grid(row=row, column=0, padx=5, pady=5, sticky=tk.W)
                ttk.Radiobutton(dialog, text="+", variable=der3_var, value=1).grid(row=row, column=1)
                ttk.Radiobutton(dialog, text="-", variable=der3_var, value=-1).grid(row=row, column=2)
                ttk.Radiobutton(dialog, text="✗", variable=der3_var, value=0).grid(row=row, column=3)
                row += 1

        def add():
            region = [xmin, xmax, monot_var.get(), convex_var.get(), der3_var.get()]
            self.regions.append(region)
        
            # Format texte selon le degré
            m = {1:"↗", -1:"↘", 0:"-"}[monot_var.get()]
            c = {1:"∪", -1:"∩", 0:"-"}[convex_var.get()] if self.degree >= 2 else "-"
            d = {1:"+", -1:"-", 0:"-"}[der3_var.get()] if self.degree >= 3 else "-"
        
            self.region_listbox.insert(tk.END, f"[{xmin:.2f}, {xmax:.2f}]  M:{m} C:{c} D3:{d}")
        
            # Visualisation
            patch = self.ax.axvspan(xmin, xmax, alpha=0.2, color='yellow')
            self.region_patches.append(patch)
            self.canvas.draw()
        
            dialog.destroy()
            self.status_var.set(f"Région ajoutée [{xmin:.2f}, {xmax:.2f}]")
    
        ttk.Button(dialog, text="Ajouter", command=add).grid(row=row, column=0, columnspan=4, pady=20)
        ttk.Button(dialog, text="Annuler", command=dialog.destroy).grid(row=row+1, column=0, columnspan=4)
        
        # ============ Construction contraintes ============
    def build_constraints_arrays(self, kn):
        """Construit les tableaux de contraintes selon le degré"""
        monot = np.zeros(kn)
        cv = np.zeros(kn + 1)
        der3 = np.zeros(kn + 1)
        if self.constraints_mode == "uniform":
        # Contraintes uniformes
            monot = np.ones(kn) * self.monot_var.get()
        
            # Pour la convexité
            if self.degree >= 2:
                cv = np.ones(kn + 1) * self.convex_var.get()
            else:
                cv = np.zeros(kn + 1)  # Pas de convexité pour degré 1
        
                # Pour la dérivée 3
                if self.degree >= 3:
                    der3 = np.ones(kn + 1) * self.der3_var.get()
                else:
                    der3 = np.zeros(kn + 1)  # Pas de dérivée 3 pour degré < 3
        
                    return monot, cv, der3
    
        else:  # Contraintes par région
           
        
            for i in range(kn):
                x1, x2 = self.knots[i], self.knots[i + 1]
            
                for r in self.regions:
                    if r[1] > x1 and r[0] < x2:
                        # Monotonie (tous degrés) - toujours appliquée
                        if r[2] != 0:
                            monot[i] = r[2]
                    
                            # Convexité - seulement si degré >= 2
                            if self.degree >= 2 and r[3] != 0:
                                cv[i] = r[3]
                                cv[i + 1] = r[3]
                    
                                # Dérivée 3 - seulement si degré >= 3
                                if self.degree >= 3 and r[4] != 0:
                                    der3[i] = r[4]
                                    der3[i + 1] = r[4]
        
        return monot, cv, der3

    # ============ Exécution ============
    def run_regression(self):
        """Lance la régression"""
        # Vérifications
        if self.xtab is None:
            messagebox.showwarning("Attention", "Chargez des données")
            return
    
        if self.knots is None:
            messagebox.showwarning("Attention", "Définissez les nœuds")
            return
    
        # Vérifier que le module approprié est disponible
        if self.degree == 1 and SplineLinearQuantile is None:
            messagebox.showerror("Erreur", "Module spline linéaire non disponible")
            return
        if self.degree == 2 and SplineQuadraticQuantile is None:
            messagebox.showerror("Erreur", "Module spline quadratique non disponible")
            return
        if self.degree == 3 and SplineCubicQuantBspkn3 is None:
            messagebox.showerror("Erreur", "Module spline cubique non disponible")
            return
        if self.degree == 4 and SplineQuantBspkn4 is None:
            messagebox.showerror("Erreur", "Module spline quartique non disponible")
            return
    
        self.progress.start()
        self.status_var.set("Calcul en cours...")
        self.root.update()
    
        try:
            # Construire les contraintes
            kn = len(self.knots) - 1
            monot, cv, der3 = self.build_constraints_arrays(kn)
        
            tau = self.tau_var.get()
            solver = self.solver_var.get()
        
            print(f"Lancement avec solveur: {solver}")
            print(f"Degré: {self.degree}")
            print(f"Monotonie: {monot}")
            if self.degree >= 2:
                print(f"Convexité: {cv}")
            if self.degree >= 3:
                    print(f"Dérivée 3: {der3}")
                    
                    # Appel à la fonction appropriée selon le degré
            if self.degree == 1:
                        self.spline_result = SplineLinearQuantile(
                            self.xtab, self.ytab, self.knots, tau, 
                            monot=monot, solver=solver
                        )
            elif self.degree == 2:
                        self.spline_result = SplineQuadraticQuantile(
                            self.xtab, self.ytab, self.knots, tau, 
                            monot=monot, cv=cv, solver=solver
                        )
            elif self.degree == 3:
                        self.spline_result = SplineCubicQuantBspkn3(
                            self.xtab, self.ytab, self.knots, tau, 
                            monot, cv, der3, solver=solver
                        )
            else:  # degree == 4
                        self.spline_result = SplineQuantBspkn4(
                            self.xtab, self.ytab, self.knots, tau, 
                            monot, cv, der3, solver=solver
                        )
        
            self.progress.stop()
        
            if self.spline_result is not None:
                            self.plot_result()
                            self.status_var.set(f"Réussi! τ={tau:.2f} (solveur: {solver})")
            else:
                            self.status_var.set("Échec - non convergence")
            
        except Exception as e:
            self.progress.stop()
            self.status_var.set(f"Erreur: {str(e)[:50]}...")
            messagebox.showerror("Erreur", str(e))
            import traceback
            traceback.print_exc()
    
    def export_results(self):
        """Exporte les résultats"""
        if self.spline_result is None:
            messagebox.showwarning("Attention", "Aucun résultat")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[
                ("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg"),
                ("CSV", "*.csv"), ("NumPy", "*.npy")
            ])
        
        if filename:
            try:
                if filename.endswith('.csv'):
                    x = np.linspace(0, 1, 500)
                    y = self.spline_result(x)
                    pd.DataFrame({'x': x, 'y': y}).to_csv(filename, index=False)
                elif filename.endswith('.npy'):
                    np.save(filename, self.spline_result.c)
                else:
                    self.fig.savefig(filename, dpi=300, bbox_inches='tight')
                self.status_var.set(f"Exporté: {filename}")
            except Exception as e:
                messagebox.showerror("Erreur", str(e))
    def generate_custom_test(self):
     """Génère des données de test avec la fonction personnalisée"""
     try:
        n = self.test_npoints_var.get()
        x = np.linspace(0, 1, n)
        
        # Évaluer la fonction
        func_str = self.test_function_var.get()
        
        # Remplacer les fonctions mathématiques
        func_str = func_str.replace('sin(', 'np.sin(')
        func_str = func_str.replace('cos(', 'np.cos(')
        func_str = func_str.replace('exp(', 'np.exp(')
        func_str = func_str.replace('log(', 'np.log(')
        func_str = func_str.replace('sqrt(', 'np.sqrt(')
        func_str = func_str.replace('pi', 'np.pi')
        func_str = func_str.replace('randn(', 'np.random.randn(')
        func_str = func_str.replace('rand(', 'np.random.rand(')
        
        # Évaluer
        y = eval(func_str, {'np': np, 'x': x,'n':n})
        
        self.xtab = x
        self.ytab = y
        self.clear_all()
        self.plot_data()
        self.data_info.set(f"{n} points (test personnalisé)")
        self.status_var.set("Données test générées")
        
     except Exception as e:
        messagebox.showerror("Erreur", f"Erreur dans la fonction: {str(e)}")

def main():
    root = tk.Tk()
    app = QuantileSplineApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
