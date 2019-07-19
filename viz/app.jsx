import React from "react";
import ReactDOM from "react-dom";
import { BrowserRouter, Route, Link } from "react-router-dom";
import Navbar from "react-bootstrap/Navbar";
import Nav from "react-bootstrap/Nav";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import { BAMViewer } from "alignment.js";

import ErrorCorrection from "./error_correction.jsx";
import SuperReadGraph from "./super_read_graph.jsx";
import QuasispeciesReconstruction from "./quasispecies_reconstruction.jsx";

import "./style.scss";


function NavLink(props) {
  return (<Nav.Link as={Link} to={props.to}>
    {props.header}
  </Nav.Link>);
}

function App() {
	return (<BrowserRouter>
		<Navbar bg="light">
			<Navbar.Brand>ACME Quasispecies Reconstruction</Navbar.Brand>
			<Nav className="mr-auto">
				<NavLink to="/" header="Mapped Reads" />
				<NavLink to="/error-correction" header="Error Correction" />
				<NavLink to="/super-read-graph" header="Super Read Graph" />
				<NavLink to="/quasispecies-reconstruction" header="Quasispecies Reconstruction" />
			</Nav>
		</Navbar>
		<Container>
			<Route
				exact
				path="/"
				component={props => <BAMViewer {...props} data_url="api/sorted.bam" />}
			/>
			<Route path="/error-correction" component={ErrorCorrection} />
			<Route path="/super-read-graph" component={SuperReadGraph} />
			<Route path="/quasispecies-reconstruction" component={QuasispeciesReconstruction} />
		</Container>
	</BrowserRouter>);
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement("div"))
);
